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
//using namespace std;

// ADOLC
// tapeless - means first derivative only - not a problem here
//#include "adolc.h"
#define ADOLC_TAPELESS // Going to use the tapeless method
//#define NUMBER_DIRECTIONS 2 // B0 and sigma2
#include "adouble.h" // only this is needed for tapeless
typedef adtl::adouble adouble; // necessary for tapeless - see manual

// Preliminaries
#ifdef WIN32
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl
#else
   #define SEXPDLLExport SEXP
#endif

// For golden search
#define resphi (2 - ((1+sqrt(5))/2))

// Function declarations
adouble fobj(adouble f, double M, double Catch, adouble bexp);
adouble fobj_grad(adouble f, double M, double Catch, adouble bexp);
adouble simpleNRtogetf(adouble f, double M, double Catch, adouble bexp);
adouble ratner_search(adouble f1, adouble f2, adouble f3, double M, double Catch, adouble bexp);

double fobj(double f, double M, double Catch, double bexp);
double ratner_search(double f1, double f2, double f3, double M, double Catch, double bexp);

void project_b(adouble* bexp, adouble* b, adouble* f, adouble** n,
                double* Catch, adouble B0,
                double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs);

void project_bh(adouble* bexp, adouble* b, adouble* h, adouble** n,
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

adouble logdnorm(double x,adouble mean,adouble sigma);

// Not actually used
adouble dnorm(double x,adouble mean,adouble sigma);

//******************************************************************************
// Golden search functions
// But these do not work when using adoubles due to conditionals
// might be easier to use a simple Newton Raphson solver using the gradient
// newf = oldf - (objfunc / grad(objfunc))
// is this going to work though? need grad fobj wrt f

// Objective function for estimating f
adouble fobj(adouble f, double M, double Catch, adouble bexp)
{
    adouble Catch_hat;
    Catch_hat = bexp*(1-exp(-f-M))*(f / (f+M));
    return ((Catch_hat - Catch)*(Catch_hat - Catch));
}

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
    for (i=0; i<100; i++)
    {
	//Rprintf("In solver f: %f\n", f.getValue());
	newf = f - fobj(f,M,Catch,bexp) / fobj_grad(f,M,Catch,bexp);
	if (sqrt(pow(f-newf,2)) > 1e-9) f = newf;
	//condassign(f,1e-9 - sqrt(pow(f-newf,2)), f,newf);
	// Need to set some kind of limit on f - pick 100
	if (f > 100) f = 100;
	//condassign(f,f-100,100); 
    }
    return f; // might need better convergence instead of just running 100 times
}
//sage: c = (f / (f+m)) * (1 - exp(-m-f)) * B
//sage: c
//-(e^(-f - m) - 1)*B*f/(f + m)
//sage: diff(c,f)
//B*f*e^(-f - m)/(f + m) - (e^(-f - m) - 1)*B/(f + m) + (e^(-f - m) - 1)*B*f/(f + m)^2

//sage: chat = (f / (f+m)) * (1 - exp(-m-f)) * B
//sage: c = var('c')
//sage: fobj = (chat - c)^2
//sage: diff(fobj,f)
//-2*((e^(-f - m) - 1)*B*f/(f + m) + c)*(B*f*e^(-f - m)/(f + m) - (e^(-f - m) - 1)*B/(f + m) + (e^(-f - m) - 1)*B*f/(f + m)^2)

// Straight double version
double fobj(double f, double M, double Catch, double bexp)
{
    double Catch_hat;
    Catch_hat = bexp*(1-exp(-f-M))*(f / (f+M));
    return ((Catch_hat - Catch)*(Catch_hat - Catch));
}

// This compiles but don't work
// conditionals?
adouble ratner_search(adouble f1, adouble f2, adouble f3, double M, double Catch, adouble bexp)
{
    // termination criteria
    if (abs(f1.getValue() - f3.getValue()) < 1e-9)
	return ((f1 +f3)/2);
    // calc new f using golden search
    adouble fnew = f2 + resphi * (f3 - f2);
    Rprintf("fnew %f\n", fnew.getValue());
    // evaluate obj func at fmid and fnew
    if (fobj(fnew,M,Catch,bexp).getValue() < fobj(f2,M,Catch,bexp).getValue())
    {
	ratner_search(f2,fnew,f3,M,Catch,bexp);
    }
    else
    {
	ratner_search(fnew,f2,f1,M,Catch,bexp);
    }
}

// Straight double version
double ratner_search(double f1, double f2, double f3, double M, double Catch, double bexp)
{
    // termination criteria
    if (abs(f1 - f3) < 1e-9)
	return ((f1 +f3)/2);
    // calc new f using golden search
    double fnew = f2 + resphi * (f3 - f2);
    //Rprintf("fnew %f\n", fnew);
    // evaluate obj func at fmid and fnew
    if (fobj(fnew,M,Catch,bexp) < fobj(f2,M,Catch,bexp))
    {
	ratner_search(f2,fnew,f3,M,Catch,bexp);
    }
    else
    {
	ratner_search(fnew,f2,f1,M,Catch,bexp);
    }
}
//******************************************************************************
// bexp
//******************************************************************************

// Project population using f estimated with Golden Search
void project_b(adouble* bexp, adouble* b, adouble* f, adouble** n,
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

  // Estimate first f
  // Double version
  //f[0] = ratner_search(0,1,100,M,Catch[0],bexp[0].getValue());
    f[0] = 1;
    f[0] = simpleNRtogetf(f[0], M, Catch[0], bexp[0]);
    Rprintf("f0: %f\n", f[0].getValue());


    // Set up SRR parameters
    alpha = (4*steepness*R0) / (5*steepness-1);
    beta = B0*(1-steepness) / (5*steepness-1);

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
	//f[yrcount] = ratner_search(0,1,100,M,Catch[yrcount],bexp[yrcount].getValue());
    f[yrcount] = simpleNRtogetf(f[yrcount-1], M, Catch[yrcount], bexp[yrcount]);
	for (i=0;i<nages;i++)
	    b[yrcount] = b[yrcount] + (n[i][yrcount] * mat[i] * wght[i]);
    }

// Print bexp and gradient
//    for (yrcount = 0; yrcount < nyrs; yrcount++)
//	Rprintf("bexp[%i} value: %f and gradient: %f \n", yrcount, bexp[yrcount].getValue(), bexp[yrcount].getADValue());

delete [] p;
}

// Project population using h as Charlie's original - not sure that this is right
void project_bh(adouble* bexp, adouble* b, adouble* h, adouble** n,
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

    // Is this safe...?
    h[0] = Catch[0] / bexp[0];
    if (h[0] < 0)
	h[0] = 0;
    // If catch > bexp, should crash out 
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

delete [] p;
}


//******************************************************************************
// Function to calculate the index hat
//******************************************************************************
void calc_index_hat(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,
                adouble* q, double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices)
{
    int i, j, nonnanyrs;
    adouble mean_log_ind_over_bexp;

    // loop over each index
    for (i=0;i<nindices;i++)
    {
	//q[i] = 0;
	mean_log_ind_over_bexp = 0;
	// Get bexp
	//project_bh(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
	project_b(bexp, b, h, n, Catch, B0, steepness, M, mat, sel, wght, amin, amax, nyrs);
	nonnanyrs = 0;
	//Rprintf("start mean log ind grad: %f\n", mean_log_ind_over_bexp.getADValue());
	//Rprintf("bexp[0] grad: %f\n", bexp[0].getADValue());
	// loop over years to then calc the mean
	for (j=0;j<nyrs;j++)
	{
	    // test for nan years in index (i.e. empty years)
	    if (!isnan(index[i][j]) & !isnan(bexp[j].getValue()))
	    {
		nonnanyrs++;
		//mean_log_ind_over_bexp = mean_log_ind_over_bexp + log(index[i][j] / bexp[j]);
		// Rewrite log section as AD getting messed up
		mean_log_ind_over_bexp = mean_log_ind_over_bexp + (log(index[i][j]) - log(bexp[j]));
		// index is fine and bexp has a grad
		// then why does the log bit fuck up the grad?
		//Rprintf("index[%i]: %f\n", j, index[i][j]);
		//Rprintf("bexp[%i]: %f\n", j, bexp[j].getValue());
		//Rprintf("bexp[%i] grad: %f\n", j, bexp[j].getADValue());
		//Rprintf("log bit with bexp[%i]: %f\n", j, (log(index[i][j] / bexp[j])).getValue());
		//Rprintf("log bit with bexp[%i] grad: %f\n", j, (log(index[i][j] / bexp[j])).getADValue());
		//Rprintf("log bit with bexp[%i] grad: %f\n", j, (log(index[i][j]) - log(bexp[j])).getADValue()); // looks OK
		//Rprintf("bit with bexp[%i] grad: %f\n", j, ((index[i][j] / bexp[j])).getADValue()); // Is OK. So maybe it's the log operation that is hurting it
		//Rprintf("Any operation on bexp[%i] grad: %f\n", j, ( log(bexp[j])).getADValue());
	    }
	}
	// get mean
	//Rprintf("almost end mean log ind grad: %f\n", mean_log_ind_over_bexp.getADValue());
	mean_log_ind_over_bexp = mean_log_ind_over_bexp / nonnanyrs;
	//Rprintf("end mean log ind grad: %f\n", mean_log_ind_over_bexp.getADValue());
	q[i] = exp(mean_log_ind_over_bexp);
	//Rprintf("q grad: %f\n", q[i].getADValue());
	//Rprintf("q: %f\n", q[i].getValue());
	// bexp and q grads look fine (compared to approx grads in R)
	// so wtf?
	for (j=0;j<nyrs;j++)
	{
	    index_hat[i][j] = bexp[j] * q[i];
//	    Rprintf("bexp[%i] grad: %f\n", j, bexp[j].getADValue());
	    //Rprintf("index_hat[%i] grad: %f\n", j, index_hat[i][j].getADValue());
	    //Rprintf("log index_hat[%i] grad: %f\n", j, (log(bexp[j]) + log(q[i])).getADValue());
//	    Rprintf("bexp [%i] * q[%i] grad: %f\n", j,j, (bexp[j]*q[i]).getADValue());
	}
    }
}

// Not actually called
adouble dnorm(double x,adouble mean,adouble sigma)
{
     return ((1 / sqrt(2*M_PI*sigma*sigma)) * (exp(-((x-mean)*(x-mean))/(2 * sigma*sigma))));
}

// Second argument is indexhat
adouble logdnorm(double x,adouble mean,adouble sigma)
{
    //Rprintf("mean grad: %f\n", mean.getADValue());
     return(log(1 / sqrt(2*M_PI*sigma*sigma)) - (((x-mean)*(x-mean))/(2 * sigma*sigma)));
}

//******************************************************************************
// Functino to calculate the logl
//******************************************************************************
adouble get_logl(adouble* bexp, adouble* b, adouble* h, adouble** n,
		adouble** index_hat, double* Catch, double** index, adouble B0,adouble sigma2,
                adouble* q, double steepness, double M, double* mat,
                double* sel, double* wght, int amin, int amax, int nyrs, int nindices)
{
    int i, j, k;
    adouble total_logl=0;
    adouble index_logl_yr=0;
    adouble temp;
//    Rprintf("start total_logl: %f\n", total_logl.getValue());
//    Rprintf("start total_logl grad: %f\n", total_logl.getADValue());
    calc_index_hat(bexp, b, h, n, index_hat, Catch, index, B0, q, steepness, M, mat, sel, wght, amin, amax, nyrs, nindices);
    // Loop over indices to get the total_logl
	//Rprintf("index_hat[0][0] grad: %f\n", index_hat[0][0].getADValue());
    for (i = 0; i<nindices; i++)
    {
	index_logl_yr = 0;
	//Rprintf("start index_logl grad: %f\n", index_logl_yr.getADValue());
	for (j=0; j<nyrs; j++)
	{
	// check for NAN in index, i.e. not all years have index data
	    if (!isnan(index[i][j]))
	    {
//		Rprintf("log index_hat[%i]: %f\n", j,log(index_hat[i][j]).getValue());
//		Rprintf("log index_hat[%i] grad: %f\n", j, log(index_hat[i][j]).getADValue());
		//Rprintf("index_hat[%i] grad: %f\n", j, (index_hat[i][j]).getADValue());
		//Rprintf("log(index / index_hat[%i]) grad: %f\n", j, log(index[i][j]/index_hat[i][j]).getADValue());
		index_logl_yr = index_logl_yr + (logdnorm(log(index[i][j]),log(index_hat[i][j]), sqrt(sigma2)));
		temp = index[i][j] / index_hat[i][j];
		Rprintf("temp %f\n", temp.getValue());
		Rprintf("temp grad %f\n", temp.getADValue());
		Rprintf("log temp %f\n", log(temp).getValue());
		Rprintf("log temp grad %f\n", log(temp).getADValue());
	    }
	}
	Rprintf("end index_logl grad: %f\n", index_logl_yr.getADValue()); // nan
	total_logl = total_logl + index_logl_yr;
    }
    Rprintf("end total_logl: %f\n", total_logl.getValue());
    Rprintf("end total_logl grad: %f\n", total_logl.getADValue()); // nan
    return total_logl;
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

//*********** Start evaluating gradients ***********************

    // Deactivate all independent adoubles
    B0.setADValue(0);
    sigma2.setADValue(0);
    // Activate B0
    B0.setADValue(1);

    total_logl = get_logl(bexp, b, f, n,
		index_hat, Catch, index, B0, sigma2,
		q, steepness, M, mat,
		sel, wght, amin, amax, nyrs, nindices);

    double B0_grad = total_logl.getADValue();

    // Deactivate all
    B0.setADValue(0);
    sigma2.setADValue(0);
    // Activate sigma2
    sigma2.setADValue(1);

    total_logl = 0;
    total_logl = get_logl(bexp, b, f, n,
		index_hat, Catch, index, B0, sigma2,
		q, steepness, M, mat,
		sel, wght, amin, amax, nyrs, nindices);

    double sigma2_grad = total_logl.getADValue();

//Rprintf("total_logl: %f\n", total_logl.getValue());
//Rprintf("sigma2_grad: %f\n", sigma2_grad);
//Rprintf("B0_grad: %f\n", B0_grad);


//******************** Return all the bits we want *********************************

    SEXP bexpSEXP, bSEXP, hSEXP, qSEXP, loglSEXP, out, outnames, loglnames, index_hatSEXP, index_dim;

    // bexp b
    PROTECT(bexpSEXP = allocVector(REALSXP,nyrs));
    PROTECT(bSEXP = allocVector(REALSXP,nyrs));
    PROTECT(hSEXP = allocVector(REALSXP,nyrs));
    for (i=0; i<nyrs; i++)
    {
	REAL(bexpSEXP)[i] = bexp[i].getValue();
	REAL(bSEXP)[i] = b[i].getValue();
	REAL(hSEXP)[i] = f[i].getValue();
    }

    // Should be one q for every index
    // q
    PROTECT(qSEXP = allocVector(REALSXP,nindices));
    for (i=0; i<nindices; i++)
	REAL(qSEXP)[i] = q[i].getValue();

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


    // Cleaning up and freeing memory
    delete [] Catch;
    delete [] mat;
    delete [] sel;
    delete [] wght;
    delete [] q;
    delete [] bexp;
    delete [] b;
    delete [] h;

    for(i=0; i<nages; i++)
	delete [] n[i];

    for(i=0; i<nindices; i++)
	delete [] index_hat[i];

    UNPROTECT(10);
    return(out);
}





