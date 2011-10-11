// FLsp
// Catch and index years can be different
// You can always get one more biomass year than catch year
// Index is compared to biomass to calculate the likelihood
// So you can have one more year in index than in catch - probably won't have but it's possible

// Moved headers and all the top gumpf to FLsp.h

//#define ADOLC_TAPELESS // Going to use the tapeless method
//#include "adouble.h" // only this is needed for tapeless
//typedef adtl::adouble adouble; // necessary for tapeless - see manual
//#include <stdlib.h>

// use Rcpp
//#include <Rcpp.h>
//#include <math.h>

#include "FLsp.h"

#define pi 3.1415926535897932384626433832795

using namespace Rcpp;

// Function headers
//void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k);
// assume error in I = q B is multiplicative and lognormal with a constant coefficient of variation
//void q_obs_mult_log(adouble* B, adouble* q, NumericVector I);
// other ways of estimating q (see Polacheck et al 1993)
// additive and normal with constant standard deviation
// additive and normal with constant coefficient of variation
//void ll_obs(adouble* B, NumericVector C, NumericVector I, adouble* q, adouble* Ihat, adouble* sigma2, adouble* ll);


RcppExport SEXP testflspCpp(SEXP C_sexp)
{
  //Rprintf("In testflspCpp\n");
  NumericVector C(C_sexp);
  C(0) = C(0) * 2;
  return wrap(C);
}

// This is the old function that is called from R using .Call
// Left in for old times sake
// Can also use the exposed functions in FLsp_exposed.cpp
RcppExport SEXP flspCpp(SEXP C_sexp, SEXP I_sexp, SEXP r_sexp, SEXP p_sexp, SEXP k_sexp)
{
  int i;

  NumericVector C(C_sexp);
  NumericVector I(I_sexp);
  //double r = as<double>(r_sexp);
  double p = as<double>(p_sexp);
  //double k = as<double>(k_sexp);

  int nyrs = C.size();
  // Make a vector for Biomass, same size as the Catch
  // Actually has one more year than catch
  //NumericVector B(clone(C));
  NumericVector B(nyrs+1);
  adouble* B_ad = new adouble[nyrs+1];
  // Index should have as many years as it wants - could be one more than catch
  // For the moment leave as same as same size - most likely
  NumericVector Ihat(nyrs);
  NumericVector res(nyrs);
  NumericVector res_grad_r(nyrs);
  NumericVector res_grad_k(nyrs);
  adouble* Ihat_ad = new adouble[nyrs];

  // We could estimate Binit, r, k and q
  // But Binit and q are impossible to estimate
  // So Binit = k and q is estimated depending on error structure
  adouble r_ad, k_ad;
  adouble* q = new adouble;
  *q = 0;
  adouble* sigma2 = new adouble;
  *sigma2 = 0;
  adouble* ll = new adouble;
  *ll = 0;

  // initialise
  r_ad = as<double>(r_sexp);
  k_ad = as<double>(k_sexp);

    double ll_grad_r;
    double ll_grad_k;

	double extinct_val = 0;

// Tapeless evaluation
// Might be quicker to do this taped rather than repeating the code
// Evaluate with r
  r_ad.setADValue(1);
  k_ad.setADValue(0);
  // initialise B0 - should be part of the solving loop?
  B_ad[0] = k_ad;
  project_biomass(B_ad, C, p, r_ad, k_ad, extinct_val);
  q_obs_mult_log(B_ad, q, I);
  ll_obs(B_ad, C, I, q, Ihat_ad, sigma2, ll);
  ll_grad_r = (*ll).getADValue();
  for (i=0; i<C.size(); i++)
    res_grad_r(i) = (I(i) - Ihat_ad[i]).getADValue();

// Evaluate with k
  r_ad.setADValue(0);
  k_ad.setADValue(1);
  // initialise B0 - should be part of the solving loop?
  B_ad[0] = k_ad;
  project_biomass(B_ad, C, p, r_ad, k_ad, extinct_val);
  q_obs_mult_log(B_ad, q, I);
  ll_obs(B_ad, C, I, q, Ihat_ad, sigma2, ll);
  ll_grad_k = (*ll).getADValue();
  for (i=0; i<C.size(); i++)
    res_grad_k(i) = (I(i) - Ihat_ad[i]).getADValue();

  for (i=0; i<C.size(); i++)
  {
    //B(i) = B_ad[i].getValue();
    Ihat(i) = Ihat_ad[i].getValue();
    res(i) = I(i) - Ihat(i);
  }

	for (i=0; i<B.size(); i++)
		B(i) = B_ad[i].getValue();
		
//return List::create(Named("B",B),
//List output.create(Named("B",B),
List output = List::create(Named("B",B),
                    Named("Ihat",Ihat),
                    Named("qhat",(*q).getValue()),
                    Named("sigma2",(*sigma2).getValue()),
                    Named("ll",(*ll).getValue()),
                    Named("ll_grad_r",ll_grad_r),
                    Named("ll_grad_k",ll_grad_k),
		    Named("res",res),
		    Named("res_grad_r",res_grad_r),
		    Named("res_grad_k",res_grad_k)
                    );

// clean up
delete[] B_ad;
delete[] Ihat_ad;
delete q;
delete sigma2;
delete ll;
                    
return output;
}

// c 10 long, 0 - 9
// b can be 11 long

void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k, double extinct_val)
{
  int yr;
	// Timing is wrong - biomass has one more year than catch
	/*
  for (yr = 1; yr<C.size(); yr++)
  {
    B[yr] = B[yr-1] + (r / p) * B[yr-1] * (1 - pow((B[yr-1] / k),p)) - C(yr-1);
		// check if population has collapsed
    B[yr] = fmax(B[yr],extinct_val);
  }
  */
  for (yr = 0; yr<C.size(); yr++)
  {
//		Rprintf("yr %i\n", yr);
//		Rprintf("C[yr] %f\n", C(yr));
//		Rprintf("B[yr] %f\n", B[yr].getValue());
    B[yr+1] = B[yr] + (r / p) * B[yr] * (1 - pow((B[yr] / k),p)) - C(yr);
		// check if population has collapsed
    B[yr+1] = fmax(B[yr+1],extinct_val);
  }
}

// Calculate qhat assuming observation error
// and assuming error in I = q B is multiplicative and lognormal with a
// constant coefficient of variation
// Iy = q By e^E    E ~ N(0, sigma2)
void q_obs_mult_log(adouble* B, adouble* q, NumericVector I)
{
  // I is same dims as C
   // Get qhat first of all
  // qhat = exp(1/n sum (log (Iy / By))
  // over years for which I has data
  int yr, n;
  *q = 0;
  n = 0;
  for (yr = 0; yr<I.size(); yr++)
    // check if index has a value
    if (!__isnan(I(yr)))
    {
      n++;
      *q = *q + log(I(yr) / B[yr]);
    }
  *q = exp(*q / n);
}

void ll_obs(adouble* B, NumericVector C, NumericVector I, adouble* q, adouble* Ihat, adouble* sigma2, adouble* ll)
{
  // q By = Iy = Cy / Ey
  // Ey is fishing effort during the year y
  
  // vhat = log(C/E)y - log(C/E)haty
  // C / E = observed catch rate
  // (C / E)hat is model predicted catch rate
  // C / E = I
  // (C / E)hat = qhat B = Ihat
  int yr, n;
  //adouble ll=0;
  adouble* vhat_ad = new adouble[C.size()];

  for (yr = 0; yr<C.size(); yr++)
    Ihat[yr] = B[yr] * *q;

  n = 0;
  *sigma2 = 0;
  for (yr = 0; yr<C.size(); yr++)
    // check if index has a value
    if (!__isnan(I(yr)))
    {
      n++;
      vhat_ad[yr] = log(I(yr) / Ihat[yr]);
      *sigma2 = *sigma2 + pow(vhat_ad[yr],2);
    }
  *sigma2 = *sigma2 / n;

  n = 0;
  *ll = 0;
  for (yr = 0; yr<C.size(); yr++)
    // check if index has a value
    if (!__isnan(I(yr)))
    {
      n++;
      *ll = *ll - pow(vhat_ad[yr],2) / (2 * *sigma2);
    }
    *ll = *ll - n * log (sqrt(2*pi* *sigma2));
    
    delete[] vhat_ad;
}

