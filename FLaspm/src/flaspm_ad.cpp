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

void pop_dyn_Edwards(adouble* bexp, adouble* b, adouble* h, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

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
  adouble alpha, beta;
  // Get the equilibrium population, initial bexp and b and alpha and beta
  pop_dyn_common(nages, nyrs, alpha, beta, n, b, bexp, M, sel, wght, mat, B0, steepness);


  // Estimate first f
  // Double version
  //f[0] = ratner_search(0,1,100,M,Catch[0],bexp[0].getValue());
    f[0] = 1;
  //  f[0] = simpleNRtogetf(f[0], M, Catch[0], bexp[0]);
    Rprintf("f0: %f\n", f[0].getValue());

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
  adouble alpha, beta;

  // Get the equilibrium population, initial bexp and b and alpha and beta
  pop_dyn_common(nages, nyrs, alpha, beta, n, b, bexp, M, sel, wght, mat, B0, steepness);


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


}


//******************************************************************************
// Entry function
//******************************************************************************
extern "C" SEXPDLLExport aspm_ad(SEXP CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
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

    for (i=0; i<nindices; i++)
    for (j=0; j<nyrs; j++)
	index[i][j] = REAL(VECTOR_ELT(indexSEXP,i))[j];

    // Set up the AD variables
    // Can just the assignment which initialises gradient to 0
    // Or can use x.setValue to just set the non-AD part
    // Or can use adouble x(1,2) to see non-AD to 1 and AD to 2
    // Or can use x.setADValue(1) to set AD value
    adouble B0 = asReal(B0SEXP);
    adouble sigma2 = asReal(sigma2SEXP);

    adouble* q = new adouble[nindices];
    //adouble* q = new adouble[1];
    adouble* bexp = new adouble[nyrs];
    adouble* b = new adouble[nyrs];
    adouble* h = new adouble[nyrs];
    adouble* f = new adouble[nyrs];

// Setting up f if not estimating - for test purposes only
//f <- c(0.67, 0.669, 0.305, 0.142, 0.091, 0.088, 0.109, 0.075)
f[0] = 0.67;
f[1] = 0.669;
f[2] = 0.305;
f[3] = 0.142;
f[4] = 0.091;
f[5] = 0.088;
f[6] = 0.109;
f[7] = 0.075;

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

    adouble total_logl = 0;


// Test pop.dyn functions
pop_dyn_Edwards(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
Rprintf("bexp[1] %f\n", bexp[1].getValue());

//******************************************************************************
// Outputs
//******************************************************************************
  return CatchSEXP;
}

