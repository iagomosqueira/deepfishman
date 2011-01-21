// Order of includes is really important
// adolc.h needs to be at or near the top else compilation breaks

#define ADOLC_TAPELESS // Going to use the tapeless method
//#define NUMBER_DIRECTIONS 2 // B0 and sigma2
#include "adouble.h" // only this is needed for tapeless
//#include "adolc.h"
typedef adtl::adouble adouble; // necessary for tapeless - see manual

// If not tapeless then use this
// #include <adolc.h> // only this seems to compile but then not tapeless

#include <stdlib.h>

//#include <ctype.h>

// R headers
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <math.h>

// Preliminaries
#ifdef WIN32
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl
#else
   #define SEXPDLLExport SEXP
#endif

// For golden search
#define resphi (2 - ((1+sqrt(5))/2))

typedef enum tagFLRConstSRR
	{
   Francis         = 1,
   Edwards         = 2
  } model;


//******************************************************************************
// Function definitions
//******************************************************************************

// Need logl for each model
adouble calc_logl_Francis(double** index, adouble* bexp, adouble* qhat,
                          int nindices, int nyrs, double indextime, double M, adouble* h);
adouble calc_logl_Edwards(double** index, adouble** index_hat, adouble sigma2,
                          int nindices, int nyrs);

// qhat is function of Bexp and index
// just qhat per index
// same function for both? Could be, but Francis uses straight mean, Edwards is geometric mean
void calc_qhat_geomean(adouble* bexp, double** index, adouble* qhat, int nindices, int nyrs,
                    double indextime, double M, adouble* f);
void calc_qhat_mean(adouble* bexp, double** index, adouble* qhat, int nindices, int nyrs,
                    double indextime, double M, adouble* f);

// index.hat
// same function for both
// index.hat = qhat * Bexp
void calc_index_hat (adouble* bexp, adouble* qhat, adouble** indexhat,
                    int nindices, int nyrs, double indextime, double M, adouble* f);

// pop.dyn
// seperate functions
// but have similar starts - getting equib values etc.
void pop_dyn_Francis(adouble* bexp, adouble* b, adouble* f, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

void pop_dyn_Edwards(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

void pop_dyn_common(int nages, int nyrs, adouble* alpha, adouble* beta, adouble** n,
                    adouble* b, adouble* bexp, double M, double* sel,
                    double* wght, double* mat, adouble B0, double steepness, double spawnlag);

adouble logdnorm(double x,adouble mean,adouble sigma);

adouble fobj(adouble f, double M, double Catch, adouble bexp);
adouble fobj_grad(adouble f, double M, double Catch, adouble bexp);
adouble simpleNRtogetf(adouble f, double M, double Catch, adouble bexp);


// logl Francis needs:
// Bmid, which comes from Bexp, which is calced by pop.dyn.Francis
// needs the index
// Calcs:
// qhat
// chat2 (sigma2)
// to get logl

// logl_Edwards()
// Estimates indices
// which calcs qhat, to calc index_hat
// to get logl

//******************************************************************************
// Common functions
//******************************************************************************
// index.hat = qhat * Bexp
void calc_index_hat (adouble* bexp, adouble* qhat, adouble** indexhat, int nindices, int nyrs,
                      double indextime, double M, adouble* f)
{
  //indextime - when does the index occur? Francis = 0.5, Edwards = 0
  int i, j;
  for (i = 0; i < nindices; i++)
    for (j =0; j < nyrs; j++)
      indexhat[i][j] = bexp[j]*exp(-indextime*(M+f[j])) * qhat[i];
}

void pop_dyn_common(int nages, int nyrs, adouble* alpha, adouble* beta, adouble** n,
                    adouble* b, adouble* bexp, double M, double* sel,
                    double* wght, double* mat, adouble B0, double steepness, double spawnlag)
{
    // Handy output for fmle
    //Rprintf("Current B0: %f\n", B0.getValue());

  // spawnlag is what point of the year the mature biomass is used for spawning
  // For Edwards it is the start of the year (0)
  // For Francis it is the end of the year (1)
  double* p = new double[nages];
  double rho = 0;
  adouble R0 = 0;
  int i;

  // Equilibrium population
  p[0] = 1;
  for (i=1; i<nages; i++)
    p[i] = p[i-1]*exp(-M);

  p[nages-1] = p[nages-1] / (1-exp(-M));

  for (i=0; i<nages; i++)
    //rho = rho + (p[i] * mat[i] * wght[i]);
    rho = rho + (p[i] * sel[i] * wght[i]);

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
  
      // Set up SRR parameters
    *alpha = (4*steepness*R0) / (5*steepness-1);
    *beta = b[0]*exp(-M*spawnlag)*(1-steepness) / (5*steepness-1);

  delete [] p;


//  Rprintf("In Common\n");
//  Rprintf("R0: %f\n", R0.getValue());
//  Rprintf("steepness: %f\n", steepness);
//  Rprintf("alpha: %f\n", (*alpha).getValue());
//  Rprintf("beta: %f\n\n", (*beta).getValue());

}

void calc_qhat_geomean(adouble* bexp, double** index, adouble* qhat, int nindices, int nyrs,
                        double indextime, double M, adouble* f)
{
//	q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
  int i, j, nonnanyrs;
  adouble mean_log_ind_over_bexp;
  for (i = 0; i<nindices; i++)
  {
    // for another index does the gradient need to be reset too?
    mean_log_ind_over_bexp = 0;
    nonnanyrs = 0;
    for (j=0;j<nyrs;j++)
    {
      // test for nan years in index (i.e. empty years)
      if (!__isnan(index[i][j]) & !__isnan(bexp[j].getValue()))
	    {
        nonnanyrs++;
		    mean_log_ind_over_bexp = mean_log_ind_over_bexp + log(index[i][j] / (bexp[j]*exp(-indextime * (M + f[j]))));
        // Rewrite log section as AD getting messed up
        //mean_log_ind_over_bexp = mean_log_ind_over_bexp + (log(index[i][j]) - log(bexp[j]));
      }
    }
  	mean_log_ind_over_bexp = mean_log_ind_over_bexp / nonnanyrs;
    //Rprintf("end mean log ind grad: %f\n", mean_log_ind_over_bexp.getADValue());
    qhat[i] = exp(mean_log_ind_over_bexp);
  }
}

void calc_qhat_mean(adouble* bexp, double** index, adouble* qhat, int nindices, int nyrs,
                    double indextime, double M, adouble* f)
{
  // surveytime - when should bexp be compared to index, i.e. when was survey done?
  // For Francis it is mid year, indextime = 0.5
  //	qhat <- apply(index[[index.count]]/bmid,c(1,6),sum,na.rm=T) / n
  int i, j, nonnanyrs;
  adouble mean_ind_over_bexp;
  for (i = 0; i<nindices; i++)
  {
    // for another index does the gradient need to be reset too?
    mean_ind_over_bexp = 0;
    nonnanyrs = 0;
    for (j=0;j<nyrs;j++)
    {
      // test for nan years in index (i.e. empty years)
      if (!__isnan(index[i][j]) & !__isnan(bexp[j].getValue()))
	    {
        nonnanyrs++;
		    mean_ind_over_bexp = mean_ind_over_bexp + (index[i][j] / (bexp[j]*exp(-indextime * (M + f[j]))));
      }
    }
    // Cast from int to double? nonnanyrs?
  	qhat[i] = mean_ind_over_bexp / nonnanyrs;
    //Rprintf("end mean log ind grad: %f\n", mean_log_ind_over_bexp.getADValue());
  }
}

// Second argument is indexhat
adouble logdnorm(double x,adouble mean,adouble sigma)
{
    //Rprintf("mean grad: %f\n", mean.getADValue());
     return(log(1 / sqrt(2*M_PI*sigma*sigma)) - (((x-mean)*(x-mean))/(2 * sigma*sigma)));
}


// Objective function for estimating f
adouble fobj(adouble f, double M, double Catch, adouble bexp)
{
    adouble Catch_hat;
    Catch_hat = bexp*(1-exp(-f-M))*(f / (f+M));
    return ((Catch_hat - Catch)*(Catch_hat - Catch));
}

// Gradient of objective function for estimating f
adouble fobj_grad(adouble f, double M, double Catch, adouble bexp)
{
    //sage_grad <- -2*((exp(-f - m) - 1)*B*f/(f + m) + c)*(B*f*exp(-f - m)/(f + m) - (exp(-f - m) - 1)*B/(f + m) + (exp(-f - m) - 1)*B*f/(f + m)^2)
    adouble grad = -2*((exp(-f - M) - 1)*bexp*f/(f + M) + Catch)*(bexp*f*exp(-f - M)/(f + M) - (exp(-f - M) - 1)*bexp/(f + M) + (exp(-f - M) - 1)*bexp*f/pow(f + M,2));
    return grad;
}

adouble simpleNRtogetf(adouble f, double M, double Catch, adouble bexp)
{
  adouble newf;
  int i;
  //Rprintf("Solving f\n");
  //Rprintf("-------------\n");
  for (i=0; i<100; i++)
  {
    newf = f - fobj(f,M,Catch,bexp) / fobj_grad(f,M,Catch,bexp);
//    Rprintf("i: %i\n", i);
//    Rprintf("newf in for loop: %f\n", newf.getValue());
//    Rprintf("newf grad in for loop: %f\n", newf.getADValue());
    if (sqrt(pow(f-newf,2)) < 1e-9) return f;
	//condassign(f,1e-9 - sqrt(pow(f-newf,2)), f,newf);
	// Need to set some kind of limit on f - pick 100
	//if (f > 100) f = 100;
	//condassign(f,f-100,100);
    f = newf;
  }
  return f; // might need better convergence instead of just running 100 times
}



//******************************************************************************
// Francis functions
//******************************************************************************

void pop_dyn_Francis(adouble* bexp, adouble* b, adouble* f, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs)
{
  // temp vars for loops and such
  int i,j, yrcount, fp;
  // more interesting vars
  int nages = amax-amin+1;
  adouble fmax = 100;
  adouble* alpha = new adouble;
  adouble* beta = new adouble;
  adouble bmat_end_last;

  // Get the equilibrium population, initial bexp and b and alpha and beta
  pop_dyn_common(nages, nyrs, alpha, beta, n, b, bexp, M, sel, wght, mat, B0, steepness, 1);


  // Estimate first f
  // Double version
  //f[0] = ratner_search(0,1,100,M,Catch[0],bexp[0].getValue());
    f[0] = 1; // Needs good first estimate of f[0] else failure to converge
    f[0] = simpleNRtogetf(f[0], M, Catch[0], bexp[0]);
    if (__isnan(f[0].getValue())) f[0] = fmax;
    //Rprintf("Initial f0 value : %f\n", f[0].getValue());
    //Rprintf("Initial f0 AD val: %f\n", f[0].getADValue() * 1e8);

// Need to check if f is something sensible
// If not set bexp to 0

    // Main population loop
    // Loop through years
  for (yrcount = 1; yrcount < nyrs; yrcount++)
  {
    // recruitment
    // 	recruitment at start of year is based on biomass at end of previous year
    // See equation for alpha in Francis paper. alpha is calced using n exp(-5M)
    // and n is already half way through year.
    bmat_end_last = b[yrcount-1]*exp(-M-f[yrcount-1]);
    n[0][yrcount] = *alpha * bmat_end_last / (*beta + bmat_end_last);
    // adult dynamics
    for (i=1; i<nages; i++)
      //n[i][yrcount] = n[i-1][yrcount-1] * exp(-M) * (1-sel[i-1] * h[yrcount-1]);
      n[i][yrcount] = n[i-1][yrcount-1] * exp(-M -(sel[i-1] * f[yrcount-1]));
    n[nages-1][yrcount] = n[nages-1][yrcount] + n[nages-1][yrcount-1] * exp(-M - (sel[nages-1] * f[yrcount-1]));
    for (i=0; i<nages; i++)
      bexp[yrcount] = bexp[yrcount] + (n[i][yrcount] * sel[i] * wght[i]);
    // Estimate f
    f[yrcount] = 1;
    //f[yrcount] = ratner_search(0,1,100,M,Catch[yrcount],bexp[yrcount].getValue());
//    Rprintf("yrcount: %i\n", yrcount);
    f[yrcount] = simpleNRtogetf(f[yrcount], M, Catch[yrcount], bexp[yrcount]);
    if (__isnan(f[yrcount].getValue())) f[yrcount] = fmax;
    for (i=0;i<nages;i++)
      b[yrcount] = b[yrcount] + (n[i][yrcount] * mat[i] * wght[i]);
  }
}

adouble calc_logl_Francis(double** index, adouble* bexp, adouble* qhat,
                          int nindices, int nyrs, double indextime, double M, adouble* h)
{
      //chat2 <- apply((index[[index.count]] / sweep(bmid,1,qhat,"*") - 1)^2,c(1,6),sum,na.rm=T) / (n-2)
	    //total.logl <- total.logl + (-n*log(sqrt(chat2)) -n*log(qhat) -apply(log(bmid[nonnaindexyears]),c(1,6),sum))
// This is the current check I am using - bad!
// if(any(bmid==0)) bmid[] <- 1e-9

  adouble chat2, total_logl, log_bmid, bmid;
  int i, j, index_count;
  total_logl = 0;
  for (i=0; i < nindices; i++)
  {
    chat2 = 0;
    index_count = 0;
    log_bmid = 0;
    for (j=0; j<nyrs; j++)
    {
      if (!__isnan(index[i][j]))
      {
        index_count++;
        bmid = bexp[j]*exp(-indextime*(M+h[j]));
        // this is dangerous
//        Rprintf("bmid: %f\n", bmid.getValue());
//        if(__isnan(bmid.getValue())) Rprintf("Yes, I am a nan\n");
        if(__isnan(bmid.getValue())) bmid = 1e-9;
        //Rprintf("bmid again: %f\n", bmid.getValue());
        //bmid = max(bmid,1e-9);
        // qhat can also be NaN
        chat2 = chat2 + pow((index[i][j] / (qhat[i] * bmid)) - 1,2);
        //Rprintf("chat2: %f", chat2.getValue());
        log_bmid = log_bmid + log(bmid);
        //Rprintf("log_bmid: %f", log_bmid.getValue());
      }
    }
    chat2 = chat2 / (index_count-2);
    total_logl = total_logl + (-index_count*log(sqrt(chat2)) - index_count*log(qhat[i]) - log_bmid);
  }
  return total_logl;
}

//******************************************************************************
// Edwards functions
//******************************************************************************
// Project population using h as Charlie's original - not sure that this is right
void pop_dyn_Edwards(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs)
{
  // temp vars for loops and such
  int i,j, yrcount, fp;
  // more interesting vars
  int nages = amax-amin+1;
  adouble* alpha = new adouble;
  adouble* beta = new adouble;
  // Get the equilibrium population, initial bexp and b and alpha and beta
  pop_dyn_common(nages, nyrs, alpha, beta, n, b, bexp, M, sel, wght, mat, B0, steepness, 0);

//Rprintf("In Edwards\n");
//Rprintf("alpha: %f\n", (*alpha).getValue());
//Rprintf("beta: %f\n", (*beta).getValue());

    // Is this safe...?
    h[0] = Catch[0] / bexp[0];
    if (h[0] < 0)
	h[0] = 0;
    // If catch > bexp, should crash out
    if (h[0] > 0.999)
	h[0] = 0.999;

    // Loop through years
    for (yrcount = 1; yrcount < nyrs; yrcount++)
    {
	// recruitment
	n[0][yrcount] = *alpha * b[yrcount-1] / (*beta + b[yrcount-1]);
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


}

adouble calc_logl_Edwards(double** index, adouble** index_hat, adouble sigma2,
                          int nindices, int nyrs)
{
  int i, j;
  adouble total_logl=0;
  adouble index_logl_yr=0;

  for (i = 0; i<nindices; i++)
  {
    index_logl_yr = 0;
    for (j=0; j<nyrs; j++)
    {
      // check for NAN in index, i.e. not all years have index data
	    if (!__isnan(index[i][j]))
	    {
        index_logl_yr = index_logl_yr + (logdnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2)));
	    }
    }
  total_logl = total_logl + index_logl_yr;
  }
  return total_logl;
}


//******************************************************************************
// Entry function
//******************************************************************************
extern "C" SEXPDLLExport aspm_ad(SEXP CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP, SEXP modelSEXP)
{

    // loop and other junk vars
    int i, j, k;

    // turn all the SEXP into something more useful
    // Do the ints
    int amin = asInteger(aminSEXP);
    int amax = asInteger(amaxSEXP);
    int nyrs = asInteger(nyrsSEXP);
    int nages = amax-amin+1;
    int model_name = asInteger(modelSEXP);

    // the scalar doubles
    double steepness = asReal(steepnessSEXP);

    // Actually, can we not just set up a pointer to the SEXP?
    // 1D vector of doubles
    // set them up
    double* Catch = new double[nyrs];
    double M = asReal(MSEXP);
    double* mat = new double[nages];
    double* sel = new double[nages];
    double* wght = new double[nages];
    // fill them in
    for (i = 0; i<nages; i++)
    {
      mat[i] = REAL(matSEXP)[i];
      sel[i] = REAL(selSEXP)[i];
     	wght[i] = REAL(wghtSEXP)[i];
    }

    for (i=0; i<nyrs; i++)
      Catch[i] = REAL(CatchSEXP)[i];

    // Set up index as 2D array
    // coming in as an FLQuants object which extends a list
    int nindices = length(indexSEXP);
    double** index = new double*[nindices];
    for (i=0;i<nindices;i++)
      index[i] = new double[nyrs];
    // And fill it in
    for (i=0; i<nindices; i++)
      for (j=0; j<nyrs; j++)
        index[i][j] = REAL(VECTOR_ELT(indexSEXP,i))[j];

    // Set up the AD variables
    // Can just the assignment which initialises gradient to 0
    // Or can use x.setValue to just set the non-AD part
    // Or can use adouble x(1,2) to see non-AD to 1 and AD to 2
    // Or can use x.setADValue(1) to set AD value

    adouble* qhat = new adouble[nindices];
    adouble* bexp = new adouble[nyrs];
    adouble* b = new adouble[nyrs];
    adouble* h = new adouble[nyrs];

    // initialise
    for (i=0;i<nyrs;i++)
    {
      bexp[i] = 0;
      b[i] = 0;
      h[i] = 0;
    }

    // Population array
    adouble** n = new adouble*[nages];
    for (i=0; i<nages; i++)
      n[i] = new adouble[nyrs];

    adouble** index_hat = new adouble*[nindices];
    for (i=0;i<nindices;i++)
      index_hat[i] = new adouble[nyrs];

    adouble total_logl = 0;
    
  double B0_grad = 0;
  double sigma2_grad = 0;


//******************************************************************************

// Model evaluations

    // Tapeless
    adouble B0 = asReal(B0SEXP);
    adouble sigma2 = asReal(sigma2SEXP);

  if (model_name == 1)
  {
    // Set up B0 deriv
    B0.setADValue(1);
    sigma2.setADValue(0);

    pop_dyn_Edwards(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
    calc_qhat_geomean(bexp, index, qhat, nindices, nyrs, 0, M, h);
    calc_index_hat (bexp, qhat, index_hat, nindices, nyrs, 0, M, h);
    total_logl = calc_logl_Edwards(index, index_hat, sigma2, nindices, nyrs);
    B0_grad = total_logl.getADValue();

    // Set up sigma2 deriv
    B0.setADValue(0);
    sigma2.setADValue(1);

    pop_dyn_Edwards(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
    calc_qhat_geomean(bexp, index, qhat, nindices, nyrs, 0, M, h);
    calc_index_hat (bexp, qhat, index_hat, nindices, nyrs, 0, M, h);
    total_logl = calc_logl_Edwards(index, index_hat, sigma2, nindices, nyrs);
    sigma2_grad = total_logl.getADValue();
  }

  //Rprintf("bexp[1] %f\n", bexp[1].getValue());
  if (model_name == 2)
  {
    // Set up B0 deriv
    B0.setADValue(1);
    sigma2.setADValue(0);
    pop_dyn_Francis(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);

//for (i=0; i<12; i++)
//{
//  Rprintf("f[%i] grad: %f\n", i, h[i].getADValue()*1e8);
//  Rprintf("bexp[%i] grad: %f\n", i, bexp[i].getADValue());
//}

    calc_qhat_mean(bexp, index, qhat, nindices, nyrs, 0.5, M, h);
    // This gradient is very wrong
    //Rprintf("qhat grad: %f\n", qhat[0].getADValue() * 1e6);
    calc_index_hat (bexp, qhat, index_hat, nindices, nyrs, 0.5, M, h);
    total_logl = calc_logl_Francis(index, bexp, qhat, nindices, nyrs, 0.5, M, h);
    B0_grad = total_logl.getADValue();
    sigma2_grad = 0;
  }


//******************************************************************************
// Outputs
//******************************************************************************
    SEXP bexpSEXP, bSEXP, hSEXP, qhatSEXP, loglSEXP, out, outnames, loglnames,
          index_hatSEXP, index_dim, nSEXP, n_dim;

    // bexp, b and h
    PROTECT(bexpSEXP = allocVector(REALSXP,nyrs));
    PROTECT(bSEXP = allocVector(REALSXP,nyrs));
    PROTECT(hSEXP = allocVector(REALSXP,nyrs));
    for (i=0; i<nyrs; i++)
    {
      REAL(bexpSEXP)[i] = bexp[i].getValue();
      REAL(bSEXP)[i] = b[i].getValue();
      REAL(hSEXP)[i] = h[i].getValue();
    }

    // qhat, should be one q for every index
    PROTECT(qhatSEXP = allocVector(REALSXP,nindices));
    for (i=0; i<nindices; i++)
      REAL(qhatSEXP)[i] = qhat[i].getValue();

    // logl - value and gradients
    PROTECT(loglSEXP = allocVector(REALSXP,3));
    REAL(loglSEXP)[0] = total_logl.getValue();
    REAL(loglSEXP)[1] = B0_grad;
    REAL(loglSEXP)[2] = sigma2_grad;

    PROTECT(loglnames = NEW_CHARACTER(3));
    SET_STRING_ELT(loglnames,0,mkChar("logl"));
    SET_STRING_ELT(loglnames,1,mkChar("logl_grad_B0"));
    SET_STRING_ELT(loglnames,2,mkChar("logl_grad_sigma2"));
    SET_NAMES(loglSEXP,loglnames);

    // index_hat
    PROTECT(index_dim     = allocVector(INTSXP, 2));
    INTEGER(index_dim)[0] = nindices;
    INTEGER(index_dim)[1] = nyrs;
    PROTECT(index_hatSEXP = Rf_allocArray(REALSXP, index_dim));
    i = 0;
    for (j = 0; j<nyrs; j++)
      for (k = 0; k <nindices; k++)
        REAL(index_hatSEXP)[i++] = index_hat[k][j].getValue();

    // n
    PROTECT(n_dim = allocVector(INTSXP, 2));
    INTEGER(n_dim)[0] = nages;
    INTEGER(n_dim)[1] = nyrs;
    PROTECT(nSEXP = Rf_allocArray(REALSXP, n_dim));
    i = 0;
    for (j = 0; j < nyrs; j++)
      for (k = 0; k < nages; k++)
        REAL(nSEXP)[i++] = n[k][j].getValue();


    // Set up the actual list to be outputted
    PROTECT(out = NEW_LIST(7));
    SET_ELEMENT(out,0,bexpSEXP);
    SET_ELEMENT(out,1,bSEXP);
    SET_ELEMENT(out,2,hSEXP);
    SET_ELEMENT(out,3,qhatSEXP);
    SET_ELEMENT(out,4,loglSEXP);
    SET_ELEMENT(out,5,index_hatSEXP);
    SET_ELEMENT(out,6,nSEXP);

    // And give it some dimnames
    PROTECT(outnames = NEW_CHARACTER(7));
    SET_STRING_ELT(outnames,0,mkChar("bexp"));
    SET_STRING_ELT(outnames,1,mkChar("bmat"));
    SET_STRING_ELT(outnames,2,mkChar("harvest"));
    SET_STRING_ELT(outnames,3,mkChar("qhat"));
    SET_STRING_ELT(outnames,4,mkChar("logl"));
    SET_STRING_ELT(outnames,5,mkChar("indexhat"));
    SET_STRING_ELT(outnames,6,mkChar("n"));
    SET_NAMES(out,outnames);


    // Cleaning up and freeing memory
    delete [] Catch;
    delete [] mat;
    delete [] sel;
    delete [] wght;
    delete [] qhat;
    delete [] bexp;
    delete [] b;
    delete [] h;

    for(i=0; i<nages; i++)
      delete [] n[i];

    for(i=0; i<nindices; i++)
      delete [] index_hat[i];

    UNPROTECT(12);
    return(out);


//  return CatchSEXP;
}

