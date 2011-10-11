// FLsp


#include <adolc.h>
#include <stdlib.h>

// use Rcpp
#include <Rcpp.h>
#include <math.h>

#define pi 3.1415926535897932384626433832795

using namespace Rcpp;

// Function headers
void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k);
// assume error in I = q B is multiplicative and lognormal with a constant coefficient of variation
void q_obs_mult_log(adouble* B, adouble* q, NumericVector I);
// other ways of estimating q (see Polacheck et al 1993)
// additive and normal with constant standard deviation
// additive and normal with constant coefficient of variation
void ll_obs(adouble* B, NumericVector C, NumericVector I, adouble* q, adouble* Ihat, adouble* sigma2, adouble* ll);


RcppExport SEXP flspCpp_tape(SEXP C_sexp, SEXP I_sexp, SEXP r_sexp, SEXP p_sexp, SEXP k_sexp)
{
  //Rprintf("In flspCpp\n");
  int i, j;

  NumericVector C(C_sexp);
  NumericVector I(I_sexp);
  //double r = as<double>(r_sexp);
  double p = as<double>(p_sexp);
  //double k = as<double>(k_sexp);

  int nyrs = C.size();
  // Make a vector for Biomass, same size as the Catch
  //NumericVector B(clone(C));
  NumericVector B(nyrs);
  NumericVector Ihat(nyrs);
  NumericVector res(nyrs);
  NumericVector res_grad_r(nyrs);
  NumericVector res_grad_k(nyrs);
    double ll_grad_r;
    double ll_grad_k;
    double ll;
    int tag=0;

  // Active variables
  adouble* B_ad = new adouble[nyrs];
  adouble* Ihat_ad = new adouble[nyrs];
  adouble r_ad, k_ad;
  adouble* rk = new adouble[2];
  adouble* q = new adouble;
  adouble* sigma2 = new adouble;
  adouble* ll_ad = new adouble;

  
  // Start the tape
    trace_on(tag);
    
  // initialise
    r_ad <<= as<double>(r_sexp);
    k_ad <<= as<double>(k_sexp);
    //rk[0] <<= as<double>(r_sexp);
    //rk[1] <<= as<double>(k_sexp);

  // initialise B0 and other stuff
  B_ad[0] = k_ad;
  *q = 0;
  *sigma2 = 0;
  *ll_ad = 0;
    project_biomass(B_ad, C, p, r_ad, k_ad);
    //project_biomass(B_ad, C, p, rk[0], rk[1]);

  q_obs_mult_log(B_ad, q, I);
  ll_obs(B_ad, C, I, q, Ihat_ad, sigma2, ll_ad);
  *ll_ad >>= ll;
  
  trace_off();

  // interrogate tape
  // get gradients (d ll / d r and d ll / dk)

  double* indeps = new double[2];
  indeps[0] = as<double>(r_sexp);
  indeps[1] = as<double>(k_sexp);
  double* grads = new double[2];

  gradient(tag,2,indeps,grads);

  // get hessian
  double** H = new double*[2];
  for (i=0; i<2; i++)
     H[i] = new double[2]; 

	hessian(tag,2,indeps,H);
  NumericMatrix Hout(2,2);
  for (i=0;i<2;i++)
      for(j=0;j<2;j++)
	  Hout(i,j) = H[i][j];

  for (i=0; i<C.size(); i++)
  {
    B(i) = B_ad[i].value();
    Ihat(i) = Ihat_ad[i].value();
    res(i) = I(i) - Ihat(i);
  }
 
//return List::create(Named("I",I));

List output = List::create(Named("B",B),
                    Named("Ihat",Ihat),
                    Named("qhat",(*q).value()),
                    Named("sigma2",(*sigma2).value()),
                    Named("ll",ll),
                    Named("ll_grad_r",grads[0]),
                    Named("ll_grad_k",grads[1]),
		    Named("hessian",Hout),
		    Named("res",res)
                    );

// clean up

// Clean pointers to adouble
delete[] B_ad;
delete[] Ihat_ad;
delete[] rk;

delete q;
delete sigma2;
delete ll_ad;

// Clean pointers to double
delete[] indeps;
delete[] grads;
for (i=0;i<2;i++)
	delete[] H[i];
delete[] H;

return output;

}


void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k)
{
  // some stuff
  //Rprintf("In project biomass\n");
  int yr;
  for (yr = 1; yr<C.size(); yr++)
  {
    B[yr] = B[yr-1] + (r / p) * B[yr-1] * (1 - pow((B[yr-1] / k),p)) - C(yr-1);
    //B[yr] = fmax(B[yr],1e-9);
    B[yr] = fmax(B[yr],0);
  }
  //Rprintf("Leaving project biomass\n");
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
  {
    // check if index has a value
    // We don't like conditionals when using tape ADOLC
    // But here we only evaluate tape once before it is written over so should be OK
    if (!__isnan(I(yr)))
    {
      n++;
      *q = *q + log(I(yr) / B[yr]);
    }
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
  {
    // check if index has a value
    if (!__isnan(I(yr)))
    {
      n++;
      vhat_ad[yr] = log(I(yr) / Ihat[yr]);
      *sigma2 = *sigma2 + pow(vhat_ad[yr],2);
    }
  }
  *sigma2 = *sigma2 / n;

  n = 0;
  *ll = 0;
  for (yr = 0; yr<C.size(); yr++)
  {
    // check if index has a value
    if (!__isnan(I(yr)))
    {
      n++;
      *ll = *ll - pow(vhat_ad[yr],2) / (2 * *sigma2);
    }
  }
    *ll = *ll - n * log (sqrt(2*pi* *sigma2));
    
    delete[] vhat_ad;
}


//******************************************************************************
// Hessian of log r and log k
RcppExport SEXP flspCpp_tape_log(SEXP C_sexp, SEXP I_sexp, SEXP r_sexp, SEXP p_sexp, SEXP k_sexp)
{
  //Rprintf("In flspCpp\n");
  int i, j;

  NumericVector C(C_sexp);
  NumericVector I(I_sexp);
  double r = as<double>(r_sexp);
  double p = as<double>(p_sexp);
  double k = as<double>(k_sexp);

  double log_r = log(r);
  double log_k = log(k);

  int nyrs = C.size();
  NumericVector B(nyrs);
  NumericVector Ihat(nyrs);
  NumericVector res(nyrs);
  NumericVector res_grad_r(nyrs);
  NumericVector res_grad_k(nyrs);
    double ll_grad_r;
    double ll_grad_k;
    double ll;
    int tag=1;

  // Active variables
  adouble* B_ad = new adouble[nyrs];
  adouble* Ihat_ad = new adouble[nyrs];
  adouble r_ad, k_ad;
  adouble log_r_ad, log_k_ad;
  adouble* q = new adouble;
  adouble* sigma2 = new adouble;
  adouble* ll_ad = new adouble;

  
  // Start the tape
    trace_on(tag);
    
  // initialise
  //  r_ad <<= as<double>(r_sexp);
  //  k_ad <<= as<double>(k_sexp);

    log_r_ad <<= log_r;
    log_k_ad <<= log_k;

    r_ad = exp(log_r_ad);
    k_ad = exp(log_k_ad);
  
  // initialise B0 and other stuff
  B_ad[0] = k_ad;
  *q = 0;
  *sigma2 = 0;
  *ll_ad = 0;
    project_biomass(B_ad, C, p, r_ad, k_ad);
    //project_biomass(B_ad, C, p, rk[0], rk[1]);

  q_obs_mult_log(B_ad, q, I);
  ll_obs(B_ad, C, I, q, Ihat_ad, sigma2, ll_ad);
  *ll_ad >>= ll;
  
  trace_off();

  // interrogate tape
  // get gradients (d ll / d r and d ll / dk)

  double* indeps = new double[2];
  indeps[0] = log_r;//as<double>(r_sexp);
  indeps[1] = log_k;//as<double>(k_sexp);
  double* grads = new double[2];

  gradient(tag,2,indeps,grads);

  //Rprintf("indeps[0] %f\n", indeps[0]);
  //Rprintf("indeps[1] %f\n", indeps[1]);
  //Rprintf("grads[0] %f\n", grads[0]);
  //Rprintf("grads[1] %f\n", grads[1]);

  // get hessian
  double** H = new double*[2];
  for (i=0; i<2; i++)
     H[i] = new double[2]; 

	hessian(tag,2,indeps,H);
  NumericMatrix Hout(2,2);
  for (i=0;i<2;i++)
      for(j=0;j<2;j++)
	  Hout(i,j) = H[i][j];

  for (i=0; i<C.size(); i++)
  {
    B(i) = B_ad[i].value();
    Ihat(i) = Ihat_ad[i].value();
    res(i) = I(i) - Ihat(i);
  }
 
//return List::create(Named("I",I));

List output = List::create(Named("B",B),
                    Named("Ihat",Ihat),
                    Named("qhat",(*q).value()),
                    Named("sigma2",(*sigma2).value()),
                    Named("ll",ll),
                    Named("ll_grad_r",grads[0]),
                    Named("ll_grad_k",grads[1]),
		    Named("hessian",Hout),
		    Named("res",res)
                    );

// clean up

// Clean pointers to adouble
delete[] B_ad;
delete[] Ihat_ad;

delete q;
delete sigma2;
delete ll_ad;

// Clean pointers to double
delete[] indeps;
delete[] grads;
for (i=0;i<2;i++)
	delete[] H[i];
delete[] H;

return output;

}




