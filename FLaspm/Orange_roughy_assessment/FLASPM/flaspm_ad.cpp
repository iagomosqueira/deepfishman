/*
Functions to calculate likelihood, gradient of likelihood,
SSB, Bexp (biomass available to fishing), Recruitment etc.

Functions from R (.Call):
  getLogl - returns logl and gradient of logl


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
#include "adouble.h" // only this is neede for tapeless
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


void index_hat(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices);

// index here has 2D, multiple indices
void logl(adouble* logl, double* Catch, double** index, adouble B0, adouble sigma2,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

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
  // right?
  adouble R0 = 0;
  adouble alpha, beta;


  // Set up equib pop
  p[0] = 1;
  for (i=1; i<nages; i++)
  {
    p[i] = p[i-1]*exp(-M);
//    Rprintf("p -1: %f\n", p[i-1]);
//    Rprintf("M -1: %f\n", M);
  }

  p[nages-1] = p[nages-1] / (1-exp(-M));
  
  for (i=0; i<nages; i++)
  {
    rho = rho + (p[i] * mat[i] * wght[i]);
//    Rprintf("p: %f\n", p[i]);
//    Rprintf("mat: %f\n", mat[i]);
//    Rprintf("wght: %f\n", wght[i]);
  }

//    Rprintf("rho: %f\n", rho);

  R0 = B0 / rho;

  for (i = 0; i<nages; i++)
  {
    n[i][0] = R0 * p[i];

  //  Rprintf("R0: %f\n", R0.getValue());
  //  Rprintf("p: %f\n", p[i]);
  //  Rprintf("n: %f\n\n", n[i][0].getValue());
  }

  // initialise
  for (i=0;i<nyrs;i++)
  {
    b[i] = 0;
    bexp[i] = 0;
//  Rprintf("bexp: %f\n", bexp[i].getValue());
//  Rprintf("bexp AD: %f\n", bexp[i].getADValue());
  }

// bexp seems OK here above

  for (i=0; i<nages; i++)
  {
    b[0] = b[0] + (n[i][0] * mat[i] * wght[i]);
    bexp[0] = bexp[0] + (n[i][0] * sel[i] * wght[i]);
//    Rprintf("sel: %f\n", sel[i]);
//    Rprintf("wght: %f\n", wght[i]);
//    Rprintf("n: %f\n\n", n[i][0].getValue());
  }
  
for (fp=0; fp<nyrs;fp++)
{
//  Rprintf("bexp: %f\n", bexp[fp].getValue());
//  Rprintf("bexp AD: %f\n", bexp[fp].getADValue());
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
      n[i][yrcount] = n[i-1][yrcount-1] * exp(-M) * (1-sel[i-1]) * h[yrcount-1];
    n[nages-1][yrcount] = n[nages-1][yrcount] + n[nages-1][yrcount-1] * exp(-M) * (1-sel[nages-1]) * h[yrcount-1];
    for (i=0; i<nages; i++)
      bexp[yrcount] = bexp[yrcount] + (n[i][yrcount] * sel[i] * wght[i]);
    h[yrcount] = Catch[yrcount] / bexp[yrcount];
  if (h[yrcount] < 0)
    h[yrcount] = 0;
  if (h[yrcount] > 0.999)
    h[yrcount] = 0.999;
  bexp[yrcount] = Catch[yrcount] / h[yrcount];

  //Rprintf("bexp: %f\n", bexp[yrcount].getValue());

  for (i=0;i<nages;i++)
    b[yrcount] = b[yrcount] + (n[i][yrcount] * mat[i] * wght[i]);
  }
  
  //double* p = new double[nages];
  delete [] p;
}


void index_hat(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices)
{

    int i, j;
//    adouble* log_ind_over_bexp = new adouble[nyrs];
//    for (j=0;j<nyrs;j++)
//	log_ind_over_bexp[j] = 0; // necessary?
    adouble mean_log_ind_over_bexp;
    adouble q = 0;

    // loop over each index
    for (i=0;i<nindices;i++)
    {
	mean_log_ind_over_bexp = 0;
	project_b(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
	//ind_over_bexp[yr] = index[i][yr] / bexp[yr] for non nan years
	// log_ind_over_bexp = log of ind_over_bexp for non nan years

	for (j=0;j<nyrs;j++)
	{
	    //log_ind_over_bexp[j] = log(index[i][j] / bexp[j]);
	    // test for nan
	    if (!isnan(index[i][j]) || !isnan(bexp[j].getValue()))
		mean_log_ind_over_bexp = mean_log_ind_over_bexp + log(index[i][j] / bexp[j]);

	}
	// mean_log_ind_over_bexp = mean(log_ind_over_bexp) for non nan years
	mean_log_ind_over_bexp = mean_log_ind_over_bexp / nyrs;
	q = exp(mean_log_ind_over_bexp);
	// q = exp(mean...)
  //q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
//  index.hat[[index.count]] <- FLQuant(q*bexp,dimnames=dm)
    }


}

//
//aspm.index <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax) {
//
//# Do index loop in this function, and return FLQuants, not in the logl function.
// for (index.count in 1:length(index))
// { 
//  bexp <- aspm.pdyn(catch,index[[index.count]],B0,hh,M,mat,sel,wght,amin,amax)
//  
//  # nuisance q
//  ind <- as.vector(index[[index.count]])
//  iyr <- as.numeric(dimnames(index[[index.count]])[['year']])
//  ys <- iyr[1]
//  yf <- iyr[length(iyr)]
//  y1 <- (ys-yr[1])+1
//  y2 <- (yf-yr[1])+1
//  q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
//  index.hat[[index.count]] <- FLQuant(q*bexp,dimnames=dm)
//  }
//  return(index.hat)
//}
//
//******************************************************************************
// Entry function
//******************************************************************************
// Could just have one entry function which returns everything
extern "C" SEXPDLLExport aspm(SEXP CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP)
{
  // loop and other junk vars
  int i, j;

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

  // return logl and gradient of log
  // and a bunch of other stuff

// set up adouble bits

// Call the functions
// test with first index
    project_b(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);

// Getting bexp back
SEXP bexpSEXP, bSEXP, hSEXP;
PROTECT(bexpSEXP = allocVector(REALSXP,nyrs));
for (i=0; i<nyrs; i++)
    REAL(bexpSEXP)[i] = bexp[i].getValue();

// gradient bits

// do some deleting and tidying up of arrays

// To get the gradient of a variable use getADValue(), e.g.
// logl.getADValue()

UNPROTECT(1);
 return bexpSEXP; 
  //return CatchSEXP;

}





/*
extern "C" SEXPDLLExport testindex(SEXP indexSEXP)
{

int i, j;
int nyrs = 11;

  int nindices = length(indexSEXP);
  Rprintf("nindices: %i\n\n", nindices);
  double** index = new double*[nindices];
  for (i=0;i<nindices;i++)
    index[i] = new double[nyrs];

  for (i=0; i<nindices; i++)
    for (j=0; j<nyrs; j++)
    {
      index[i][j] = REAL(VECTOR_ELT(indexSEXP,i))[j]; // Not sure about this one...
      Rprintf("index[i][j]: %f\n", index[i][j]);
    }

return indexSEXP;

}
*/
