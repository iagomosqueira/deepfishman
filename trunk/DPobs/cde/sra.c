
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
void pop_dyn(double,double*,double*,double*);
void pop_index(double,double*,double*,double*,double*);
double pop_index_sigma(double*,double*);

// global variables
int ymin,ymax,amin,amax,nyr,nag;
double hh,*H,*M;
double *mat,*wght,*sel;
double *C,*I;
double alp,bet;

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
  nyr  = ymax - ymin + 2;
  
  ///////////////////
  // PRELIMINARIES //
  ///////////////////
 
  // catch
  C = REAL(R_catch);

  // index
  I = REAL(R_index);

  // steepness
  hh = NUMERIC_VALUE(R_hh);
 
  // natural mortality
  M = REAL(R_M);

  // maturity-at-age
  mat = REAL(R_mat);

  // weight-at-age
  wght = REAL(R_wght);

  // selectivity
  sel = REAL(R_sel);
  
  ///////////////
  // RUN MODEL //
  ///////////////

  // doubles
  double B0;
  double sigma,sigma2;
  double nLogLk;

  double *B     = new double[nyr];  
  double *Bexp  = new double[nyr];  
  double *H     = new double[nyr];    
  double *Ipred = new double[nyr];

  // initialise doubles
  B0 = NUMERIC_VALUE(R_B0);
  
  // run model
  // :: initialise likelihood
  nLogLk = 0.;
  
  // :: predicted index
  pop_index(B0,B,Bexp,H,Ipred);
  
  // :: calculate sigma
  sigma = pop_index_sigma(I,Ipred);

  // :: calculate negative log-likelihood
  sigma2 = pow(sigma,2);
  for(y=0;y<nyr;y++) {
    if(Ipred[y]>0. && I[y]>0.) {
      nLogLk += log(sigma2) + pow(I[y]-Ipred[y],2.)/sigma2;
    }
  }

  ////////////
  // Output //
  ////////////
  int p;
  int iAge,iYear;
  SEXP dag,dyr,names,dimnames,params,srnames;

  // create dimension names
  PROTECT(dag = allocVector(INTSXP,nag));
  for(iAge=amin,a=0;a<nag;iAge++,a++) {
    INTEGER(dag)[a] = iAge;
  }

  PROTECT(dyr = allocVector(INTSXP,nyr));
  for(iYear=ymin,y=0;y<nyr;iYear++,y++) {
    INTEGER(dyr)[y] = iYear;
  }

  // B
  SEXP _B;
  PROTECT(_B = allocVector(REALSXP,nyr));

  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_B)[p++] = B[y];
  }

  setAttrib(_B,R_NamesSymbol,dyr);

  // Bexp
  SEXP _Bexp;
  PROTECT(_Bexp = allocVector(REALSXP,nyr));

  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_Bexp)[p++] = Bexp[y];
  }

  setAttrib(_Bexp,R_NamesSymbol,dyr);

  // H
  SEXP _H;
  PROTECT(_H = allocVector(REALSXP,nyr));
 
  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_H)[p++] = H[y];
  }

  setAttrib(_H,R_NamesSymbol,dyr);

  // Ipred
  SEXP _Ipred;
  PROTECT(_Ipred = allocVector(REALSXP,nyr));
 
  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_Ipred)[p++] = Ipred[y];
  }

  setAttrib(_Ipred,R_NamesSymbol,dyr);

  // sigma
  SEXP _sigma;
  PROTECT(_sigma = allocVector(REALSXP,1));
  REAL(_sigma)[0] = sigma;

  // -logLk
  SEXP _nLogLk;
  PROTECT(_nLogLk = allocVector(REALSXP,1));
  REAL(_nLogLk)[0] = nLogLk;
  
  // SR par
  SEXP _srpar;
  PROTECT(_srpar = allocVector(REALSXP,2));
  REAL(_srpar)[0] = alp;
  REAL(_srpar)[1] = bet;
  PROTECT(srnames = allocVector(STRSXP,2));
  SET_STRING_ELT(srnames,0,mkChar("alpha"));
  SET_STRING_ELT(srnames,1,mkChar("beta"));
  setAttrib(_srpar,R_NamesSymbol,srnames);

  // create combined output list
  SEXP out;
  PROTECT(out = allocVector(VECSXP,7));
  SET_VECTOR_ELT(out,0,_B);
  SET_VECTOR_ELT(out,1,_Bexp);
  SET_VECTOR_ELT(out,2,_H);
  SET_VECTOR_ELT(out,3,_Ipred);
  SET_VECTOR_ELT(out,4,_sigma);
  SET_VECTOR_ELT(out,5,_nLogLk);
  SET_VECTOR_ELT(out,6,_srpar);

  // assign names
  PROTECT(names = allocVector(STRSXP,7));
  SET_STRING_ELT(names,0,mkChar("B"));
  SET_STRING_ELT(names,1,mkChar("Bexp"));
  SET_STRING_ELT(names,2,mkChar("H"));
  SET_STRING_ELT(names,3,mkChar("Ipred"));
  SET_STRING_ELT(names,4,mkChar("sigma"));
  SET_STRING_ELT(names,5,mkChar("nLogLk"));
  SET_STRING_ELT(names,6,mkChar("srpar"));
  setAttrib(out,R_NamesSymbol,names);

  UNPROTECT(22);

  // clean up
  delete[] B,Bexp,H,Ipred;

  UNPROTECT(2);
  return out;
  
}

///////////////
// FUNCTIONS //
///////////////

void pop_dyn(double B0,double *B,double *Bexp,double *H) { 

  int a,y;

  double rho = 0.;
  double *P = new double[nag];

  double R0;
  double **N = new double*[nag];
  for(a=0;a<nag;a++) {
    N[a] = new double[nyr];
  }

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
  
  // Loop through the years

  for(y=1;y<nyr;y++) {
  
    // biomass
    B[y-1] = 0.;
    Bexp[y-1] = 0.;
    for(a=0;a<nag;a++) {
      B[y-1] += N[a][y-1] * mat[a] * wght[a];
      Bexp[y-1] += N[a][y-1] * sel[a] * wght[a] * exp(-M[a]/2.);
    }

    // harvest
    H[y-1] = 1.;
    if(Bexp[y-1]>0.) H[y-1] = C[y-1] / Bexp[y-1];
    if(H[y-1]<0.)    H[y-1] = 0.;
    if(H[y-1]>1.)    H[y-1] = 1.;
    if(H[y-1]>0.)    Bexp[y-1] = C[y-1] / H[y-1];

    // recruitment

    N[0][y] = alp * B[y-1]/(bet + B[y-1]);

    // adult dynamics

    for(a=1;a<nag;a++)
      N[a][y] = N[a-1][y-1]*exp(-M[a-1])*(1-sel[a-1]*H[y-1]);
    N[nag-1][y] = N[nag-1][y] + N[nag-1][y-1]*exp(-M[nag-1])*(1-sel[nag-1]*H[y-1]);

  }
  
  // current biomass
  B[nyr-1] = 0.;
  Bexp[nyr-1] = 0.;
  for(a=0;a<nag;a++) {
    B[nyr-1] += N[a][nyr-1] * mat[a] * wght[a];
    Bexp[nyr-1] += N[a][nyr-1] * sel[a] * wght[a] * exp(-M[a]/2.);
  }
  
  // current H
  H[nyr-1] = 0.;
  
  // clean up
  delete[] P;

  for(a=0;a<nag;a++) {
    delete[] N[a];
  }
  delete[] N;

}

void pop_index(double B0,double *B,double *Bexp,double *H,double *Ipred) {

  int y;
  int n = 0;

  double q = 0.;
  double tmp = 0.;
  
  pop_dyn(B0,B,Bexp,H);

  for(y=0;y<nyr;y++) {
    if(I[y]>0. && Bexp[y]>0.) {
      tmp += log(I[y]/Bexp[y]);
      n++;
    }
  }

  if(n>0) q = exp(tmp/n);

  for(y=0;y<nyr;y++) {
    Ipred[y] = q*Bexp[y];
  }
}

double pop_index_sigma(double *I,double *Ipred) {

  int y;
  int n = 0;

  double tmp = 0.;
  
  for(y=0;y<nyr;y++) {
    if(I[y]>0. && Ipred[y]>0. && !ISNA(Ipred[y])) {
      tmp += pow(log(I[y]/Ipred[y]),2);
      n++;
    }
  }

  if(n>0) { return(sqrt(tmp/n));
  } else  { return(0.);
  }

}



