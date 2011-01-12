// Order of includes is really important
// adolc.h needs to be at or near the top else compilation breaks

#define ADOLC_TAPELESS // Going to use the tapeless method
//#define NUMBER_DIRECTIONS 2 // B0 and sigma2
#include "adouble.h" // only this is needed for tapeless
typedef adtl::adouble adouble; // necessary for tapeless - see manual

// If not tapeless then use this
// #include <adolc.h> // only this seems to compile but then not tapeless

#include <math.h>
#include <stdlib.h>
#include <ctype.h>

// R headers
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

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
adouble calc_logl_Francis();
adouble calc_logl_Edwards();

// qhat is function of Bexp and index
// just qhat per index
// same function for both? Could be, but Francis uses straight mean, Edwards is geometric mean
void calc_qhat_geomean(adouble* bexp, double** index, adouble qhat);
void calc_qhat_mean(adouble* bexp, double** index, adouble qhat);

// index.hat
// same function for both
// index.hat = qhat * Bexp
void calc_indexhat (adouble* bexp, adouble qhat, adouble* indexhat);

// pop.dyn
// seperate functions
// but have similar starts - getting equib values etc.
void pop_dyn_Francis(adouble* bexp, adouble* b, adouble* f, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

void pop_dyn_Edwards();

void pop_dyn_common(int nages, int nyrs, adouble alpha, adouble beta, adouble** n,
                    adouble* b, adouble* bexp, double M, double* sel,
                    double* wght, double* mat, adouble B0, double steepness);


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

void pop_dyn_common(int nages, int nyrs, adouble alpha, adouble beta, adouble** n,
                    adouble* b, adouble* bexp, double M, double* sel,
                    double* wght, double* mat, adouble B0, double steepness)
{
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
    alpha = (4*steepness*R0) / (5*steepness-1);
    beta = B0*(1-steepness) / (5*steepness-1);

  delete [] p;

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
//  double* p = new double[nages];
//  double rho = 0;
  // As R0 is calced using B0, it also has to be adouble
//  adouble R0 = 0;
  adouble alpha, beta;

  // Set up equib pop
//  p[0] = 1;
//  for (i=1; i<nages; i++)
//    p[i] = p[i-1]*exp(-M);

//  p[nages-1] = p[nages-1] / (1-exp(-M));

//  for (i=0; i<nages; i++)
    //rho = rho + (p[i] * mat[i] * wght[i]);
//    rho = rho + (p[i] * sel[i] * wght[i]);

//  R0 = B0 / rho;

  // Initial population vector in year 0
//  for (i = 0; i<nages; i++)
//    n[i][0] = R0 * p[i];

  // initialise bexp and b
//  for (i=0;i<nyrs;i++)
//  {
//    b[i] = 0;
//    bexp[i] = 0;
//  }

//  for (i=0; i<nages; i++)
//  {
//    b[0] = b[0] + (n[i][0] * mat[i] * wght[i]);
//    bexp[0] = bexp[0] + (n[i][0] * sel[i] * wght[i]);
//  }

  // Get the equilibrium population, initial bexp and b and alpha and beta
  pop_dyn_common(nages, nyrs, alpha, beta, n, b, bexp, M, sel, wght, mat, B0, steepness);


  // Estimate first f
  // Double version
  //f[0] = ratner_search(0,1,100,M,Catch[0],bexp[0].getValue());
    f[0] = 1;
  //  f[0] = simpleNRtogetf(f[0], M, Catch[0], bexp[0]);
    Rprintf("f0: %f\n", f[0].getValue());


    // Set up SRR parameters
//    alpha = (4*steepness*R0) / (5*steepness-1);
//    beta = B0*(1-steepness) / (5*steepness-1);

    // Main population loop
    // Loop through years
    for (yrcount = 1; yrcount < nyrs; yrcount++)
    {
	// recruitment
	n[0][yrcount] = alpha * b[yrcount-1] / (beta + b[yrcount-1]);
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
  //f[yrcount] = simpleNRtogetf(f[yrcount-1], M, Catch[yrcount], bexp[yrcount]);
	for (i=0;i<nages;i++)
	    b[yrcount] = b[yrcount] + (n[i][yrcount] * mat[i] * wght[i]);
    }

// Print bexp and gradient
//    for (yrcount = 0; yrcount < nyrs; yrcount++)
//	Rprintf("bexp[%i} value: %f and gradient: %f \n", yrcount, bexp[yrcount].getValue(), bexp[yrcount].getADValue());

//delete [] p;

}


//******************************************************************************
// Entry function
//******************************************************************************
extern "C" SEXPDLLExport aspm_ad(SEXP CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP)
{

  return CatchSEXP;
}

