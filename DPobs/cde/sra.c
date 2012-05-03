
// ASPM
// ::
// Charles T.T. Edwards, Imperial College
// Feb 2010

#include <iostream>
#include <fstream>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;

// global functions
void pop_dyn(double*,double,double,double*,double*,double*,double*,double*);
double q_calc(double*,double*);
void index_calc(double*,double,double*);
double sigma_calc(double*,double*);

// global variables
int ymin,ymax,amin,amax,nyr,nag,lpen;

extern "C" SEXP fit(SEXP R_B0,SEXP R_catch,SEXP R_index,SEXP R_hh,SEXP R_M,SEXP R_mat,SEXP R_sel,SEXP R_wght,SEXP R_amin,SEXP R_amax,SEXP R_ymin,SEXP R_ymax) {

  //local variables
  int a,y;

  //////////////
  // Coercion //
  //////////////

  PROTECT(R_B0     = AS_NUMERIC(R_B0));
  PROTECT(R_catch  = AS_NUMERIC(R_catch));
  PROTECT(R_index  = AS_NUMERIC(R_index));
  PROTECT(R_hh     = AS_NUMERIC(R_hh));
  PROTECT(R_M      = AS_NUMERIC(R_M));
  PROTECT(R_mat    = AS_NUMERIC(R_mat));
  PROTECT(R_sel    = AS_NUMERIC(R_sel));
  PROTECT(R_wght   = AS_NUMERIC(R_wght));
  PROTECT(R_amin   = AS_NUMERIC(R_amin));
  PROTECT(R_amax   = AS_NUMERIC(R_amax));
  PROTECT(R_ymin   = AS_NUMERIC(R_ymin));
  PROTECT(R_ymax   = AS_NUMERIC(R_ymax));
  
  ////////////////
  // DIMENSIONS //
  ////////////////

  amin = INTEGER_VALUE(R_amin);
  amax = INTEGER_VALUE(R_amax);
  nag  = amax - amin + 1;

  ymin = INTEGER_VALUE(R_ymin);
  ymax = INTEGER_VALUE(R_ymax);
  nyr  = ymax - ymin + 1;
  
  ///////////////////
  // PRELIMINARIES //
  ///////////////////
  
  // initialise SSB
  double B0 = NUMERIC_VALUE(R_B0);
  
  // catch
  double* C = REAL(R_catch);
  
  // index
  double* I = REAL(R_index);

  // steepness
  double hh = NUMERIC_VALUE(R_hh);
 
  // natural mortality
  double* M = REAL(R_M);

  // maturity-at-age
  double* mat = REAL(R_mat);

  // weight-at-age
  double* wght = REAL(R_wght);

  // selectivity
  double* sel = REAL(R_sel);
  
  UNPROTECT(12);
  
  ///////////////
  // RUN MODEL //
  ///////////////
  
  // run model
  lpen = 0.;
  double* Bexp = new double[nyr];
  pop_dyn(Bexp,B0,hh,M,mat,wght,sel,C);
  //Rprintf("Bexp: ");
  //for(y=0;y<nyr;y++) Rprintf("%f; ",Bexp[y]);
  //Rprintf("\n");
  // :: catchability
  double q = q_calc(I,Bexp);
  //Rprintf("q: %f\n",q);
  // :: predicted index
  double* Ipred = new double[nyr];
  index_calc(Ipred,q,Bexp);
  //Rprintf("Ipred: ");
  //for(y=0;y<nyr;y++) Rprintf("%f; ",Ipred[y]);
  //Rprintf("\n");
  // :: calculate sigma
  double sigma = sigma_calc(I,Ipred);
  //Rprintf("sigma: %f\n",q);
  // :: calculate negative log-likelihood
  double nLogLk = 0.;
  for(y=0;y<nyr;y++) {
    if(Ipred[y]>0. && I[y]>0. && !ISNA(Ipred[y]) && !ISNA(I[y])) {
      nLogLk += log(sigma) + 0.5 * pow(log(I[y]/Ipred[y])/sigma,2.);
      //Rprintf("%i: %f; %f; %f\n",y,I[y],Ipred[y],nLogLk);
    }
  }
  // :: add penalty
  nLogLk += lpen;
  
  //Rprintf("lengths\n C:%i;\n I:%i;\n M:%i;\n mat:%i;\n wght:%i;\n sel:%i;\n\n nag:%i;\n nyr:%i\n",LENGTH(R_catch),LENGTH(R_index),LENGTH(R_M),LENGTH(R_mat),LENGTH(R_wght),LENGTH(R_sel),nag,nyr);
  
  ////////////
  // Output //
  ////////////

  // -logLk
  SEXP _nLogLk;
  PROTECT(_nLogLk = allocVector(REALSXP,1));
  REAL(_nLogLk)[0] = nLogLk;

  // clean up
  delete[] mat;
  delete[] wght;
  delete[] sel;
  delete[] M;
  delete[] C;
  delete[] I;
  delete[] Ipred;
  delete[] Bexp;
  
  UNPROTECT(1);
  //return out;
  return _nLogLk;
  
}

///////////////
// FUNCTIONS //
///////////////

void pop_dyn(double* Bexp,double B0,double hh,double* M,double* mat,double* wght,double* sel,double* C) { 

  int a,y;
  
  double alp,bet;

  double rho = 0.;
  double *P = new double[nag];

  double R0;
  double **N = new double*[nag];
  for(a=0;a<nag;a++) {
    N[a] = new double[nyr];
  }
  
  double* B    = new double[nyr];
  double* H    = new double[nyr];

  // set up eqm population
  P[0]=1;
  for(a=1;a<nag;a++)
    P[a] = P[a-1]*exp(-M[a-1]);
  P[nag-1] = P[nag-1]/(1-exp(-M[nag-1]));
  for(a=0;a<nag;a++)
    rho += P[a] * mat[a] * wght[a];
  R0 = B0 / rho;
  for(a=0;a<nag;a++)
    N[a][0] = R0 * P[a];
  
  // set up S-R parameters

  alp = (4*hh*R0)/(5*hh-1);
  bet = B0*(1-hh)/(5*hh-1);
  
  // Initial biomass
  B[0] = 0.;
  Bexp[0] = 0.;
  for(a=0;a<nag;a++) {
    B[0] += N[a][0] * mat[a] * wght[a];
    Bexp[0] += N[a][0] * sel[a] * wght[a] * exp(-M[a]/2.);
  }

  // Initial harvest
  H[0] = 1.;
  if(Bexp[0]>0.) H[0]    = C[0] / Bexp[0];
  if(H[0]<0.)    H[0]    = 0.;
  if(H[0]>1.)    H[0]    = 1.;
  if(H[0]>0.)    Bexp[0] = C[0] / H[0];

  
  // Loop through the years
  for(y=1;y<nyr;y++) {
  
    // recruitment
    N[0][y] = alp * B[y-1]/(bet + B[y-1]);

    // adult dynamics
    for(a=1;a<nag;a++)
      N[a][y] = N[a-1][y-1]*exp(-M[a-1])*(1-sel[a-1]*H[y-1]);
    N[nag-1][y] = N[nag-1][y] + N[nag-1][y-1]*exp(-M[nag-1])*(1-sel[nag-1]*H[y-1]);
    
    // biomass
    B[y] = 0.;
    Bexp[y] = 0.;
    for(a=0;a<nag;a++) {
      B[y] += N[a][y] * mat[a] * wght[a];
      Bexp[y] += N[a][y] * sel[a] * wght[a] * exp(-M[a]/2.);
    }
    
    // harvest
    H[y] = 1.;
    if(Bexp[y]>0.) H[y]    = C[y] / Bexp[y];
    if(H[y]<0.)    H[y]    = 0.;
    if(H[y]>1.)    H[y]    = 1.;
    if(H[y]>0.)    Bexp[y] = C[y] / H[y];
    
    // ll penalty
    if(H[y]==1.) lpen += 100.;
  }

  // clean up
  delete[] P;

  for(a=0;a<nag;a++) {
    delete[] N[a];
  }
  delete[] N;
  delete[] H;
  delete[] B;

}

void index_calc(double* Ipred,double q,double* Bexp) {

  int y;

  for(y=0;y<nyr;y++) {
    Ipred[y] = q*Bexp[y];
  }
}

double q_calc(double* I,double* Bexp) {

  int y;
  int n = 0;

  double tmp = 0.;

  for(y=0;y<nyr;y++) {
    if(I[y]>0. && Bexp[y]>0. && !ISNA(I[y])) {
      tmp += log(I[y]/Bexp[y]);
      n++;
    }
  }

  if(n>0) { return(exp(tmp/n));
  } else {  return(0.);
  }

}

double sigma_calc(double *I,double *Ipred) {

  int y;
  int n = 0;

  double tmp = 0.;
  
  for(y=0;y<nyr;y++) {
    if(I[y]>0. && Ipred[y]>0. && !ISNA(I[y]) && !ISNA(Ipred[y])) {
      tmp += pow(log(I[y]/Ipred[y]),2);
      n++;
    }
  }

  if(n>0) { return(sqrt(tmp/n));
  } else  { return(0.);
  }

}



