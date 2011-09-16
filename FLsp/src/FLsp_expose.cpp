// Can we expose a function?
//#include <Rcpp.h>
//using namespace Rcpp; // necessary!

#include "FLsp.h"

// Testing basic functions
double finorm(double x, double y){
	return sqrt(x*x + y*y);
}

NumericVector expose_vector(NumericVector C)
{
	return(C);
}

double expose_vector2(NumericVector A, NumericVector B)
{
	int i;
	double C=0;
	for (i=0; i< A.size(); i++)
		C = C + A(i) + B(i);
	return(C);
}

// Can we call the Project function?
//void project_biomass(adouble* B, NumericVector C, double p, adouble r, adouble k);
NumericVector project_biomassC(NumericVector C, double r, double k, double p, double extinct_val)
{
  int nyrs = C.size();
  int i;
  // Make a vector for Biomass, same size as the Catch
  NumericVector B(nyrs);
  adouble* B_ad = new adouble[nyrs];
	adouble r_ad, k_ad;
    // initialise
  r_ad = r;
  k_ad = k;
  B_ad[0] = k_ad;
  project_biomass(B_ad, C, p, r_ad, k_ad, extinct_val);
  for (i=0; i<nyrs; i++)
		B(i) = B_ad[i].getValue();

	delete[] B_ad;

	return B;
}

// A function that gets all the bits together
// This can then be called by other C functions to just return the bits we want?
List eval_FLsp(NumericVector C, NumericVector I, double r, double k, double p, double extinct_val)
{
  int i;
  int nyrs = C.size();
  double ll_grad_r;
  double ll_grad_k;

  // Make a vector for Biomass, same size as the Catch
  //NumericVector B(clone(C));
  NumericVector B(nyrs);
  adouble* B_ad = new adouble[nyrs];
  NumericVector Ihat(nyrs);
  NumericVector res(nyrs);
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
  r_ad = r;
  k_ad = k;

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
  {
    B(i) = B_ad[i].getValue();
    Ihat(i) = Ihat_ad[i].getValue();
    res(i) = I(i) - Ihat(i);
  }

	NumericVector ll_grads(2);
	ll_grads(0) = ll_grad_r;
	ll_grads(1) = ll_grad_k;

	List output = List::create(Named("B",B),
                    Named("Ihat",Ihat),
                    Named("qhat",(*q).getValue()),
                    Named("sigma2",(*sigma2).getValue()),
                    Named("ll",(*ll).getValue()),
//                    Named("ll_grad_r",ll_grad_r),
//                    Named("ll_grad_k",ll_grad_k),
										Named("ll_grads", ll_grads),
								    Named("res",res)
                    );

	// clean up
	delete[] B_ad;
	delete[] Ihat_ad;
	delete q;
	delete sigma2;
	delete ll;

	return output;
}

// Just a function to get the likelihood
//double get_logl(NumericVector C, NumericVector I, double p, double r, double k)
//{
//	List all = eval_FLsp(C, I, p, r, k);
//	double ll = as<double>(all["ll"]);
//	return ll;
//}

// Just another function to get the likelihood
// Note different argument list to match optim call
// No overloading?
double get_logl(NumericVector params, NumericVector C, NumericVector I, double p, double extinct_val)
{
	// k is second params
	// r is first
	List all = eval_FLsp(C, I, params(0), params(1), p, extinct_val);
	double ll = as<double>(all["ll"]);
	return ll;
}


NumericVector get_loglgrads(NumericVector params, NumericVector C, NumericVector I, double p, double extinct_val)
{
	List all = eval_FLsp(C, I, params(0), params(1), p, extinct_val);
	return all["ll_grads"];
}



RCPP_MODULE(yada){
	function( "finorm", &finorm);
	function( "expose_vector", &expose_vector);
	function( "expose_vector2", &expose_vector2);
	//function( "project_biomass", &project_biomassC);
	function( "project_biomass", &project_biomassC, List::create( _["C"], _["r"], _["k"], _["p"]=1.0, _["extinct_val"]=0));
	function( "eval_FLsp", &eval_FLsp, List::create( _["C"], _["I"], _["r"], _["k"],_["p"]=1.0, _["extinct_val"]=0));
//	function( "get_logl", &get_logl, List::create( _["C"], _["I"], _["p"]=1.0, _["r"], _["k"]));
	function( "get_logl", &get_logl, List::create(_["params"], _["C"], _["I"], _["p"]=1.0,_["extinct_val"]=0));
	function( "get_loglgrads", &get_loglgrads, List::create( _["params"], _["C"], _["I"], _["p"]=1.0, _["extinct_val"]=0));
}



// Need to add something to Namespace?
// exportPattern("^[[:alpha:]]+")

