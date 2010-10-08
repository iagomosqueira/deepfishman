/*
Functions to calculate likelihood, gradient of likelihood,
SSB, Bexp (biomass available to fishing), Recruitment etc.

Functions from R (.Call):

*/

// Includes
// Standard ones
#include <stdlib.h>
#include <math.h>

// R headers
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h> //?
using namespace std;

// ADOLC
// tapeless - means first derivative only - not a problem here
//#include "adolc.h"
#define ADOLC_TAPELESS // Going to use the tapeless method
//#define NUMBER_DIRECTIONS 2 // B0 and sigma2
#include "adouble.h" // only this is needed for tapeless
typedef adtl::adouble adouble; // necessary for tapeless - see manual

// junk?
//#include <iostream>
//#include <fstream>
//#include <R.h>

using namespace std;

// Preliminaries
#ifdef WIN32
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl
#else
   #define SEXPDLLExport SEXP
#endif

// Function declarations

void project_b(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);


void calc_index_hat(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,
                adouble* q,double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices);

// index here has 2D, multiple indices
adouble get_logl(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,adouble sigma2,
                adouble* q, double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices);

adouble dnorm(double x,adouble mean,adouble sigma);
adouble logdnorm(double x,adouble mean,adouble sigma);

//******************************************************************************
// bexp
//******************************************************************************

void project_b(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs)
{
  // temp vars for loops and such
  int i,j, yrcount, fp;
  // more interesting vars
  int nages = amax-amin+1;
  double* p = new double[nages];
  double rho = 0;
  // As R0 is calced using B0, it also has to be adouble
  adouble R0 = 0;
  adouble alpha, beta;

  // Set up equib pop
  p[0] = 1;
  for (i=1; i<nages; i++)
    p[i] = p[i-1]*exp(-M);

  p[nages-1] = p[nages-1] / (1-exp(-M));
  
  for (i=0; i<nages; i++)
    rho = rho + (p[i] * mat[i] * wght[i]);

  R0 = B0 / rho;

    // Initial population vector in year 0
  for (i = 0; i<nages; i++)
    n[i][0] = R0 * p[i];

    // initialise bexp and b
  for (i=0;i<nyrs;i++)
  {
    b[i] = 0;
    bexp[i] = 0;
  }

  for (i=0; i<nages; i++)
  {
    b[0] = b[0] + (n[i][0] * mat[i] * wght[i]);
    bexp[0] = bexp[0] + (n[i][0] * sel[i] * wght[i]);
  }
 
  h[0] = Catch[0] / bexp[0];
  if (h[0] < 0)
    h[0] = 0;
  if (h[0] > 0.999)
    h[0] = 0.999;

  // Set up SRR parameters
  alpha = (4*steepness*R0) / (5*steepness-1);
  beta = B0*(1-steepness) / (5*steepness-1);

  // Loop through years
  for (yrcount = 1; yrcount < nyrs; yrcount++)
  {
    // recruitment
    n[0][yrcount] = alpha * b[yrcount-1] / (beta + b[yrcount-1]);
    // adult dynamics
    for (i=1; i<nages; i++)
      n[i][yrcount] = n[i-1][yrcount-1] * exp(-M) * (1-sel[i-1] * h[yrcount-1]);
    n[nages-1][yrcount] = n[nages-1][yrcount] + n[nages-1][yrcount-1] * exp(-M) * (1-sel[nages-1] * h[yrcount-1]);
    for (i=0; i<nages; i++)
      bexp[yrcount] = bexp[yrcount] + (n[i][yrcount] * sel[i] * wght[i]);
    h[yrcount] = Catch[yrcount] / bexp[yrcount];
    if (h[yrcount] < 0)
    h[yrcount] = 0;
    if (h[yrcount] > 0.999)
    h[yrcount] = 0.999;
    bexp[yrcount] = Catch[yrcount] / h[yrcount];

  for (i=0;i<nages;i++)
    b[yrcount] = b[yrcount] + (n[i][yrcount] * mat[i] * wght[i]);
  }

  //double* p = new double[nages];
  delete [] p;
}


void calc_index_hat(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,
                adouble* q, double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices)
{

    int i, j, nonnanyrs;
//    adouble* log_ind_over_bexp = new adouble[nyrs];
//    for (j=0;j<nyrs;j++)
//	log_ind_over_bexp[j] = 0; // necessary?
    adouble mean_log_ind_over_bexp;
    //adouble q = 0;
    q[0] = 0;

    // loop over each index
    for (i=0;i<nindices;i++)
    {
	mean_log_ind_over_bexp = 0;
	project_b(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
	//ind_over_bexp[yr] = index[i][yr] / bexp[yr] for non nan years
	// log_ind_over_bexp = log of ind_over_bexp for non nan years

	    nonnanyrs = 0;
	for (j=0;j<nyrs;j++)
	{
	    //log_ind_over_bexp[j] = log(index[i][j] / bexp[j]);
	    // test for nan
	    if (!isnan(index[i][j]) & !isnan(bexp[j].getValue()))
	    {
		nonnanyrs++;
		mean_log_ind_over_bexp = mean_log_ind_over_bexp + log(index[i][j] / bexp[j]);
	    }

	}
	// mean_log_ind_over_bexp = mean(log_ind_over_bexp) for non nan years
	mean_log_ind_over_bexp = mean_log_ind_over_bexp / nonnanyrs;
	//Rprintf("mean_log_ind_over_bexp: %f\n", mean_log_ind_over_bexp.getValue());
	q[0] = exp(mean_log_ind_over_bexp);
	//Rprintf("q: %f\n", q.getValue());
	for (j=0;j<nyrs;j++)
	    index_hat[i][j] = bexp[j] * q[0];
	// q = exp(mean...)
	//for (j=0;j<nyrs;j++)
	//    Rprintf("index_hat: %f\n", index_hat[i][j].getValue()); // looks OK

	//for (j=0;j<nyrs;j++)
	//    Rprintf("bexp: %f\n", bexp[j].getValue()); // looks OK
	//bexp looks OK,
	//q and index_hat are different
  //q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
//  index.hat[[index.count]] <- FLQuant(q*bexp,dimnames=dm)
    }
}

//sum(dnorm(log(index[[index.count]]),log.indexhat.quants[[index.count]], sqrt(sigma2), TRUE),na.rm=T)
// dnorm(x,mean,sd)
// dn = (1 / (sqrt(2*pi*sigma^2))) * (exp(-(x-mean)^2/(2*sigma^2)))

adouble dnorm(double x,adouble mean,adouble sigma)
{
    //Rprintf("1/sqrt() %f\n", (1 / sqrt(2*M_PI*sigma*sigma)).getValue());
    //Rprintf("(x-mean)^2/2sig^2 %f\n",(-((x-mean)*(x-mean))/(2 * sigma*sigma)).getValue() );
    //Rprintf("out: %f\n\n",((1 / sqrt(2*M_PI*sigma*sigma)) * (exp(-((x-mean)*(x-mean))/(2 * sigma*sigma)))).getValue() );
     return ((1 / sqrt(2*M_PI*sigma*sigma)) * (exp(-((x-mean)*(x-mean))/(2 * sigma*sigma))));
}

adouble logdnorm(double x,adouble mean,adouble sigma)
{
    //Rprintf("1/sqrt() %f\n", (1 / sqrt(2*M_PI*sigma*sigma)).getValue());
    //Rprintf("(x-mean)^2/2sig^2 %f\n",(-((x-mean)*(x-mean))/(2 * sigma*sigma)).getValue() );
    //Rprintf("out: %f\n\n",((1 / sqrt(2*M_PI*sigma*sigma)) * (exp(-((x-mean)*(x-mean))/(2 * sigma*sigma)))).getValue() );
     return(log(1 / sqrt(2*M_PI*sigma*sigma)) - (((x-mean)*(x-mean))/(2 * sigma*sigma)));
}

adouble get_logl(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,adouble sigma2,
                adouble* q, double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices)
{
    int i, j, k;
//  Rprintf("In get_logl B0: %f\n", B0.getValue());
//  Rprintf("In get_logl B0 AD: %f\n", B0.getADValue());
    adouble total_logl=0;
    adouble index_logl_yr=0;
    calc_index_hat(bexp, b, h, n, index_hat, Catch, index, B0, q, steepness, M, mat, sel, wght, amin, amax, nyrs, nindices);

    // log the index_hat quants
//for (i = 0; i<nindices; i++)
 //   for(j=0; j<nyrs; j++)
    // log_index_hat[i][j] = log(index_hat[i][j])
// Or just log as we go along - otherwise we have to make a new var

for (i = 0; i<nindices; i++)
{
	index_logl_yr = 0;
    for (j=0; j<nyrs; j++)
    {
//	Rprintf("sigma2: %f\n", sigma2.getValue());
//	Rprintf("index: %f\n", index[i][j]);
//	Rprintf("index_hat: %f\n", index_hat[i][j].getValue());
	//Rprintf("dnorm bit: %f\n\n",log(dnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2))).getValue() );
	//Rprintf("dnorm bit: %f\n\n",(logdnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2))).getValue() );


	// check for isnan again
	if (!isnan(index[i][j]))
	    //index_logl_yr = index_logl_yr + log(dnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2)));
	    index_logl_yr = index_logl_yr + (logdnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2)));
    }
    // index_logl = 
	//sum(dnorm(log(index[[index.count]]),log.indexhat.quants[[index.count]], sqrt(sigma2), TRUE),na.rm=T)
    total_logl = total_logl + index_logl_yr;
// dnorm(x,mean,sd)
// dn = (1 / (sqrt(2*pi*sigma^2))) * (exp(-(x-mean)^2/(2*sigma^2)))
//adouble dnorm(double x,adouble mean,adouble sigma)
}
//Rprintf("total_logl: %f\n", total_logl.getValue());
return total_logl;
}
/*
  logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {

    # Get the FLQuants object with the estimated indices
    indexhat.quants <- aspm.index(catch,index,B0,hh,M,mat,sel,wght,amin,amax)

    log.indexhat.quants <- lapply(indexhat.quants,function(index) window(log(index),start=dims(index)$minyear,end=dims(index)$maxyear))
    total.logl <- 0
    for (index.count in 1:length(index))
	total.logl <- total.logl + sum(dnorm(log(index[[index.count]]),log.indexhat.quants[[index.count]], sqrt(sigma2), TRUE),na.rm=T)
    return(total.logl)
  }
  */
//******************************************************************************
// Entry function
//******************************************************************************
// Could just have one entry function which returns everything
extern "C" SEXPDLLExport aspm(SEXP CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP)
{
  // loop and other junk vars
  int i, j, k;

  // turn all the SEXP into something more useful
  // Do the ints
  int amin = asInteger(aminSEXP);
  int amax = asInteger(amaxSEXP);
  int nyrs = asInteger(nyrsSEXP);
  int nages = amax-amin+1;
  
  // the scalar doubles
  double steepness = asReal(steepnessSEXP);

  // Actually, can we not just set up a pointer to the SEXP?
  // 1D vector of doubles
  // set them up
  double* Catch = new double[nyrs];
  //double* M = new double[nages];
    double M = asReal(MSEXP);
  double* mat = new double[nages];
  double* sel = new double[nages];
  double* wght = new double[nages];
  // fill them in
  for (i = 0; i<nages; i++)
  {
   // M[i] = REAL(MSEXP)[i];
    mat[i] = REAL(matSEXP)[i];
    sel[i] = REAL(selSEXP)[i];
    wght[i] = REAL(wghtSEXP)[i];
  }

  for (i=0; i<nyrs; i++)
  {
    Catch[i] = REAL(CatchSEXP)[i];
  }

  // Set up index as 2D array
  // coming in as an FLQuants object which extends a list
  int nindices = length(indexSEXP);
  double** index = new double*[nindices];
  for (i=0;i<nindices;i++)
    index[i] = new double[nyrs];

  for (i=0; i<nindices; i++)
    for (j=0; j<nyrs; j++)
      index[i][j] = REAL(VECTOR_ELT(indexSEXP,i))[j]; // Not sure about this one...

  // Set up the AD variables
  // Can just the assignment which initialises gradient to 0
  // Or can use x.setValue to just set the non-AD part
  // Or can use adouble x(1,2) to see non-AD to 1 and AD to 2
  // Or can use x.setADValue(1) to set AD value
  adouble B0 = asReal(B0SEXP);
  adouble sigma2 = asReal(sigma2SEXP);

  adouble* q = new adouble[1];
  adouble* bexp = new adouble[nyrs];
  adouble* b = new adouble[nyrs];
  adouble* h = new adouble[nyrs];

  for (i=0;i<nyrs;i++)
  {
    bexp[i] = 0;
    b[i] = 0;
    h[i] = 0;
  }
  adouble** n = new adouble*[nages];
  for (i=0; i<nages; i++)
    n[i] = new adouble[nyrs];

  adouble** index_hat = new adouble*[nindices];
  for (i=0;i<nindices;i++)
    index_hat[i] = new adouble[nyrs];

  //for (i=0; i<nindices; i++)
  //  for (j=0; j<nyrs; j++)
  //    index_hat[i][j] = 0;
  
  // return logl and gradient of log
  // and a bunch of other stuff

// set up adouble bits

// Call the functions
// test with first index
//project_b(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
//
//calc_index_hat(bexp, b, h, n,
//		index_hat, Catch, index, B0,
//                q, steepness, M, mat,
//                sel, wght, amin, amax, nyrs, nindices);
adouble total_logl = 0;

// Deactivate all
  B0.setADValue(0);
  sigma2.setADValue(0);
// Activate B0
// what is value of B0
  B0.setADValue(1);
//  Rprintf("B0: %f\n", B0.getValue());
//  Rprintf("B0 AD: %f\n", B0.getADValue());


total_logl = 0;
total_logl = get_logl(bexp, b, h, n,
		index_hat, Catch, index, B0, sigma2,
                q, steepness, M, mat,
                sel, wght, amin, amax, nyrs, nindices);

double B0_grad = total_logl.getADValue();
//Rprintf("total_logl: %f\n", total_logl.getValue());
//Rprintf("total_logl AD: %f\n", total_logl.getADValue());

// Deactivate all
  B0.setADValue(0);
  sigma2.setADValue(0);
// Activate sigma2
  sigma2.setADValue(1);

total_logl = 0;
total_logl = get_logl(bexp, b, h, n,
		index_hat, Catch, index, B0, sigma2,
                q, steepness, M, mat,
                sel, wght, amin, amax, nyrs, nindices);

double sigma2_grad = total_logl.getADValue();
//Rprintf("total_logl: %f\n", total_logl.getValue());
//Rprintf("total_logl AD: %f\n", total_logl.getADValue());

// OUPUTS
// Need to return:
// bexp
// b
// h
// q
// logl
// logl grad
//
// Make a list

SEXP bexpSEXP, bSEXP, hSEXP, qSEXP, loglSEXP, out, outnames, loglnames;
// bexp b
PROTECT(bexpSEXP = allocVector(REALSXP,nyrs));
PROTECT(bSEXP = allocVector(REALSXP,nyrs));
PROTECT(hSEXP = allocVector(REALSXP,nyrs));
for (i=0; i<nyrs; i++)
{
    REAL(bexpSEXP)[i] = bexp[i].getValue();
    REAL(bSEXP)[i] = b[i].getValue();
    REAL(hSEXP)[i] = h[i].getValue();

}

PROTECT(qSEXP = allocVector(REALSXP,1));
REAL(qSEXP)[0] = q[0].getValue();
PROTECT(loglSEXP = allocVector(REALSXP,3));
REAL(loglSEXP)[0] = total_logl.getValue();
//REAL(loglSEXP)[1] = total_logl.getADValue(0);
//REAL(loglSEXP)[2] = total_logl.getADValue(1); // the other grad
REAL(loglSEXP)[1] = B0_grad;
REAL(loglSEXP)[2] = sigma2_grad;

PROTECT(loglnames = NEW_CHARACTER(3));
SET_STRING_ELT(loglnames,0,mkChar("logl"));
SET_STRING_ELT(loglnames,1,mkChar("logl_grad_B0"));
SET_STRING_ELT(loglnames,2,mkChar("logl_grad_sigma2"));
SET_NAMES(loglSEXP,loglnames);

SEXP index_hatSEXP;
SEXP index_dim;
PROTECT(index_dim     = allocVector(INTSXP, 2));       
INTEGER(index_dim)[0] = nindices;
INTEGER(index_dim)[1] = nyrs;

PROTECT(index_hatSEXP = Rf_allocArray(REALSXP, index_dim)); 
i = 0;
for (j = 0; j<nyrs; j++)
    for (k = 0; k <nindices; k++)
	REAL(index_hatSEXP)[i++] = index_hat[k][j].getValue();


// Set up the actual list to be outputted
PROTECT(out = NEW_LIST(6));
SET_ELEMENT(out,0,bexpSEXP);
SET_ELEMENT(out,1,bSEXP);
SET_ELEMENT(out,2,hSEXP);
SET_ELEMENT(out,3,qSEXP);
SET_ELEMENT(out,4,loglSEXP);
SET_ELEMENT(out,5,index_hatSEXP);

// And give it some dimnames
PROTECT(outnames = NEW_CHARACTER(6));
SET_STRING_ELT(outnames,0,mkChar("Bexp"));
SET_STRING_ELT(outnames,1,mkChar("B"));
SET_STRING_ELT(outnames,2,mkChar("h"));
SET_STRING_ELT(outnames,3,mkChar("q"));
SET_STRING_ELT(outnames,4,mkChar("logl"));
SET_STRING_ELT(outnames,5,mkChar("indexhat"));
SET_NAMES(out,outnames);

// gradient bits

// do some deleting and tidying up of arrays

// To get the gradient of a variable use getADValue(), e.g.
// logl.getADValue()

UNPROTECT(10);
return(out);
//return index_hatSEXP;
// return bexpSEXP; 
  //return CatchSEXP;

}





