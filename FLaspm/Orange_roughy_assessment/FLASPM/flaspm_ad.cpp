/*
Functions to calculate likelihood, gradient of likelihood,
SSB, Bexp (biomass available to fishing), Recruitment etc.

Functions from R (.Call):
  getLogl - returns logl and gradient of logl


*/

// Includes
// Also need to include ADOLC
// tapeless?
//#define ADOLC_TAPELESS
//#include "c:/adolc-1.10.1/adolc/adouble.h"
//typedef adtl::adouble adouble;

#include "adolc.h"

#include <stdlib.h>
#include <math.h>

using namespace std;



//#include <iostream>
//#include <fstream>
//#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h> //?

using namespace std;

// Preliminaries
#ifdef WIN32
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl
#else
   #define SEXPDLLExport SEXP
#endif

// Function declarations

void bexp(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, double* index, adouble B0,
                double steepness, double* M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);


void index_hat(adouble* index_hat, double* Catch, double* index, adouble B0,
                double steepness, double* M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

// index here has 2D, multiple indices
void logl(adouble* logl, double* Catch, double** index, adouble B0, adouble sigma2,
                double steepness, double* M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);



void bexp(adouble* bexp, double* Catch, double* index, adouble B0,
                double steepness, double* M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

//******************************************************************************
// bexp
//******************************************************************************

void bexp(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, double* index, adouble B0,
                double steepness, double* M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs)
{
  // temp vars for loops and such
  int i,j, yrcount;
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
    p[i] = p[i-1]*exp(-M[i-1]);
  p[nages-1] = p[nages-1] / (1-exp(-M[nages-1]));
  
  for (i=0; i<nages; i++)
    rho = rho + (p[i] * mat[i] * wght[i]);

  R0 = B0 / rho;

  for (i = 0; i<nages; i++)
    n[i][0] = R0 * p[i];


  // initialise
  for (i=0;i<nages;i++)
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
      n[i][yrcount] = n[i-1][yrcount-1] * exp(-M[i]) * (1-sel[i-1]) * h[yrcount-1];
    n[nages-1][yrcount] = n[nages-1][yrcount] + n[nages-1][yrcount-1] * exp(-M[nages-1]) * (1-sel[nages-1]) * h[yrcount-1];
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
  double* Catch = new double[nages];
  double* M = new double[nages];
  double* mat = new double[nages];
  double* sel = new double[nages];
  double* wght = new double[nages];
  // fill them in
  for (i = 0; i<nages; i++)
  {
    Catch[i] = REAL(CatchSEXP)[i];
    M[i] = REAL(MSEXP)[i];
    mat[i] = REAL(matSEXP)[i];
    sel[i] = REAL(selSEXP)[i];
    wght[i] = REAL(wghtSEXP)[i];
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

// Other variables we need
  adouble* bexp = new adouble[nyrs];
  adouble* b = new adouble[nyrs];
  adouble* h = new adouble[nyrs];
// Initialise these properly!
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

// gradient bits

// do some deleting and tidying up of arrays





  
  return CatchSEXP;

}





/*
extern "C" SEXPDLLExport SepVPA_ad(SEXP xStock, SEXP xControl, SEXP xRefHarvest)
   {
   //Input VPA Control
   FLashVPA VPA(xStock);

   //Run Seperable VPA
   VPA.SepVPA(xControl, xRefHarvest);

   //Display Ns and Fs
   return VPA.Return();
   }
*/

