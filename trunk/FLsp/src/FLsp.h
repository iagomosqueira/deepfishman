#define ADOLC_TAPELESS // Going to use the tapeless method
#include "adouble.h" // only this is needed for tapeless
typedef adtl::adouble adouble; // necessary for tapeless - see manual
#include <stdlib.h>

// use Rcpp
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k, double extinct_val);
// assume error in I = q B is multiplicative and lognormal with a constant coefficient of variation
void q_obs_mult_log(adouble* B, adouble* q, NumericVector I);
// other ways of estimating q (see Polacheck et al 1993)
// additive and normal with constant standard deviation
// additive and normal with constant coefficient of variation
void ll_obs(adouble* B, NumericVector C, NumericVector I, adouble* q, adouble* Ihat, adouble* sigma2, adouble* ll);




