
// ASPM
// ::
// Charles T.T. Edwards, Imperial College
// Feb 2010

#define ADOLC_TAPELESS
#include "c:/adolc-1.10.1/adolc/adouble.h"
typedef adtl::adouble adouble;

#include <iostream>
#include <fstream>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;

// global functions
void pop_dyn(adouble,adouble*,adouble*,adouble*);
void pop_index(adouble,adouble*,adouble*,adouble*,adouble*);

// global variables
int ymin,ymax,ymin_obj,ymax_obj,amin,amax,nyr,nyr_obj,nag;
double hh,*F,*M;
double *mat,*wght,*sel;
double *Cb,*I;

extern "C" SEXP fit(SEXP R_B0,SEXP R_sigma2,SEXP R_catch,SEXP R_index,SEXP R_hh,SEXP R_M,SEXP R_mat,SEXP R_sel,SEXP R_wght,SEXP R_amin,SEXP R_amax,SEXP R_ymin,SEXP R_ymax,SEXP R_ymin_obj,SEXP R_ymax_obj) {

  //local variables
  int a,y;

  //////////////
  // Coercion //
  //////////////

  PROTECT(R_B0     = AS_NUMERIC(R_B0));
  PROTECT(R_sigma2 = AS_NUMERIC(R_sigma2));
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
  PROTECT(R_ymin_obj   = AS_NUMERIC(R_ymin_obj));
  PROTECT(R_ymax_obj   = AS_NUMERIC(R_ymax_obj));
  
  ////////////////
  // DIMENSIONS //
  ////////////////

  amin = INTEGER_VALUE(R_amin);
  amax = INTEGER_VALUE(R_amax);
  nag  = amax - amin + 1;

  ymin = INTEGER_VALUE(R_ymin);
  ymax = INTEGER_VALUE(R_ymax);
  nyr  = ymax - ymin + 1;

  ymin_obj = INTEGER_VALUE(R_ymin_obj);
  ymax_obj = INTEGER_VALUE(R_ymax_obj);
  nyr_obj  = ymax_obj - ymin_obj + 1;
  
  ///////////////////
  // PRELIMINARIES //
  ///////////////////
 
  // catch
  Cb = REAL(R_catch);

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

  // adoubles
  adouble B0_ad;
  adouble sigma2_ad;
  adouble nLogLk_ad;

  adouble *B_ad     = new adouble[nyr];  
  adouble *Bexp_ad  = new adouble[nyr];  
  adouble *F_ad     = new adouble[nyr];    
  adouble *Ipred_ad = new adouble[nyr_obj];

  // output vector
  SEXP Out = R_NilValue;
  PROTECT(Out = NEW_NUMERIC(3));

  // initialise adoubles
  B0_ad     = NUMERIC_VALUE(R_B0);
  sigma2_ad = NUMERIC_VALUE(R_sigma2);
  
  // Tapeless loop - evaluate the function with each individual AD variable 'active'

  // 'deactivate' the independent variables
  B0_ad.setADValue(0);
  sigma2_ad.setADValue(0);

  // 'activate' the appropriate independent variable
  B0_ad.setADValue(1); 

  // run model
  // :: initialise likelihood
  nLogLk_ad = 0;
  
  // :: predicted index
  pop_index(B0_ad,B_ad,Bexp_ad,F_ad,Ipred_ad);

  // :: calculate negative likelihood
  for(y=0;y<nyr_obj;y++) {
    if(Ipred_ad[y]>0. && I[y]>0.) {
      nLogLk_ad -= log(sigma2_ad) + pow(I[y]-Ipred_ad[y],2)/sigma2_ad;
    }
  }
  //Rprintf("nLogLk: %f; ",nLogLk_ad.getValue());
  //Rprintf("gradient: %f\n",nLogLk_ad.getADValue());

  // unload AD dependent variable (gradient)
  REAL(Out)[0] = nLogLk_ad.getADValue();
  
  // 'deactivate' the independent variables
  B0_ad.setADValue(0);
  sigma2_ad.setADValue(0);

  // 'activate' the appropriate independent variable
  sigma2_ad.setADValue(1); 

  // run model
  // :: initialise likelihood
  nLogLk_ad = 0;

  // :: predicted index
  pop_index(B0_ad,B_ad,Bexp_ad,F_ad,Ipred_ad);
    
  // :: calculate negative likelihood
  for(y=0;y<nyr_obj;y++) {
    if(Ipred_ad[y]>0. && I[y]>0.) {
      nLogLk_ad -= log(sigma2_ad) + pow(I[y]-Ipred_ad[y],2)/sigma2_ad;
    }
  }
  //Rprintf("nLogLk: %f; ",nLogLk_ad.getValue());
  //Rprintf("gradient: %f\n",nLogLk_ad.getADValue());

  // unload AD dependent variable (gradient)
  REAL(Out)[1] = nLogLk_ad.getADValue();

  // get likelihood
  REAL(Out)[2] = nLogLk_ad.getValue();

  //Rprintf("nLogLk: %f: ",nLogLk_ad.getValue());
  //Rprintf("\n");

  ////////////
  // Output //
  ////////////
  int p;
  int iAge,iYear;
  SEXP dag,dyr,dyr_obj,names,dimnames,params;

  // create dimension names
  PROTECT(dag = allocVector(INTSXP,nag));
  for(iAge=amin,a=0;a<nag;iAge++,a++) {
    INTEGER(dag)[a] = iAge;
  }

  PROTECT(dyr = allocVector(INTSXP,nyr));
  for(iYear=ymin,y=0;y<nyr;iYear++,y++) {
    INTEGER(dyr)[y] = iYear;
  }

  PROTECT(dyr_obj = allocVector(INTSXP,nyr_obj));
  for(iYear=ymin_obj,y=0;y<nyr_obj;iYear++,y++) {
    INTEGER(dyr_obj)[y] = iYear;
  }

  // B
  SEXP _B;
  PROTECT(_B = allocVector(REALSXP,nyr));

  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_B)[p++] = B_ad[y].getValue();
  }

  setAttrib(_B,R_NamesSymbol,dyr);

  // Bexp
  SEXP _Bexp;
  PROTECT(_Bexp = allocVector(REALSXP,nyr));

  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_Bexp)[p++] = Bexp_ad[y].getValue();
  }

  setAttrib(_Bexp,R_NamesSymbol,dyr);

  // F
  SEXP _F;
  PROTECT(_F = allocVector(REALSXP,nyr));
 
  p=0;
  for(y=0;y<nyr;y++) {
    REAL(_F)[p++] = F_ad[y].getValue();
  }

  setAttrib(_F,R_NamesSymbol,dyr);

  // Ipred
  SEXP _Ipred;
  PROTECT(_Ipred = allocVector(REALSXP,nyr_obj));
 
  p=0;
  for(y=0;y<nyr_obj;y++) {
    REAL(_Ipred)[p++] = Ipred_ad[y].getValue();
  }

  setAttrib(_Ipred,R_NamesSymbol,dyr_obj);

  // grad
  SEXP _grad;
  PROTECT(_grad = allocVector(REALSXP,2));
  REAL(_grad)[0] = REAL(Out)[0];
  REAL(_grad)[1] = REAL(Out)[1];
  PROTECT(params = allocVector(STRSXP,2));
  SET_STRING_ELT(params,0,mkChar("B0"));
  SET_STRING_ELT(params,1,mkChar("sigma2"));
  setAttrib(_grad,R_NamesSymbol,params);

  // -logLk
  SEXP _nLogLk;
  PROTECT(_nLogLk = allocVector(REALSXP,1));
  REAL(_nLogLk)[0] = REAL(Out)[2];

  // create combined output list
  SEXP out;
  PROTECT(out = allocVector(VECSXP,6));
  SET_VECTOR_ELT(out,0,_B);
  SET_VECTOR_ELT(out,1,_Bexp);
  SET_VECTOR_ELT(out,2,_F);
  SET_VECTOR_ELT(out,3,_Ipred);
  SET_VECTOR_ELT(out,4,_grad);
  SET_VECTOR_ELT(out,5,_nLogLk);

  // assign names
  PROTECT(names = allocVector(STRSXP,6));
  SET_STRING_ELT(names,0,mkChar("B"));
  SET_STRING_ELT(names,1,mkChar("Bexp"));
  SET_STRING_ELT(names,2,mkChar("F"));
  SET_STRING_ELT(names,3,mkChar("Ipred"));
  SET_STRING_ELT(names,4,mkChar("grad"));
  SET_STRING_ELT(names,5,mkChar("nLogLk"));
  setAttrib(out,R_NamesSymbol,names);

  UNPROTECT(27);

  // clean up
  delete[] B_ad,Bexp_ad,F_ad,Ipred_ad;

  UNPROTECT(1);
  return out;
  
}

///////////////
// FUNCTIONS //
///////////////

void pop_dyn(adouble B0_ad,adouble *B_ad,adouble *Bexp_ad,adouble *F_ad) { 

  int a,y;

  double rho = 0.;
  double *P = new double[nag];

  adouble R0_ad;
  adouble **C_ad = new adouble*[nag];
  adouble **N_ad = new adouble*[nag];
  for(a=0;a<nag;a++) {
    C_ad[a] = new adouble[nyr];
    N_ad[a] = new adouble[nyr];
  }

  // set up eqm population
  P[0]=1;
  for(a=1;a<nag;a++)
    P[a] = P[a-1]*exp(-M[a]);
  P[nag-1] = P[nag-1]/(1-exp(-M[nag-1]));
  for(a=0;a<nag;a++)
    rho += P[a] * mat[a] * wght[a];
  R0_ad = B0_ad / rho;
  for(a=0;a<nag;a++)
    N_ad[a][0] = R0_ad * P[a];

  B_ad[0] = 0.;
  Bexp_ad[0] = 0.;
  for(a=0;a<nag;a++) {
    B_ad[0] += N_ad[a][0] * mat[a] * wght[a];
    Bexp_ad[0] += N_ad[a][0] * sel[a] * wght[a];
  }

  F_ad[0] = Cb[0] / Bexp_ad[0];
  if(F_ad[0]<0.) F_ad[0] = 0.;
  if(F_ad[0]>1.) F_ad[0] = 0.999;
  Bexp_ad[0] = Cb[0] / F_ad[0];
  
  // set up S-R parameters

  adouble alp_ad = (4*hh*R0_ad)/(5*hh-1);
  adouble bet_ad = B0_ad*(1-hh)/(5*hh-1);
  
  // Loop through the years

  for(y=1;y<nyr;y++) {

    // recruitment

    N_ad[0][y] = alp_ad * B_ad[y-1]/(bet_ad + B_ad[y-1]);

    // adult dynamics

    for(a=1;a<nag;a++)
      N_ad[a][y] = N_ad[a-1][y-1]*exp(-M[a])*(1-sel[a-1]*F_ad[y-1]);
    N_ad[nag-1][y] = N_ad[nag-1][y] + N_ad[nag-1][y-1]*exp(-M[nag-1])*(1-sel[nag-1]*F_ad[y-1]);

    B_ad[y] = 0.;
    Bexp_ad[y] = 0.;
    for(a=0;a<nag;a++) {
      B_ad[y] += N_ad[a][y] * mat[a] * wght[a];
      Bexp_ad[y] += N_ad[a][y] * sel[a] * wght[a];
    }

    F_ad[y] = Cb[y] / Bexp_ad[y];
    if(F_ad[y]<0.) F_ad[y] = 0.;
    if(F_ad[y]>1.) F_ad[y] = 0.999;
    Bexp_ad[y] = Cb[y] / F_ad[y];
  }

  // clean up
  delete[] P;

  for(a=0;a<nag;a++) {
    delete[] C_ad[a];
    delete[] N_ad[a];
  }
  delete[] C_ad;
  delete[] N_ad;

}

void pop_index(adouble B0_ad,adouble *B_ad,adouble *Bexp_ad,adouble *F_ad,adouble *Ipred_ad) {

  int y,y_obj;
  int n = 0;

  adouble q_ad = 0.;
  adouble tmp_ad = 0.;
  
  pop_dyn(B0_ad,B_ad,Bexp_ad,F_ad);

  for(y_obj=0,y=ymin_obj-ymin;y_obj<nyr_obj;y_obj++,y++) {
    if(I[y_obj]>0. && Bexp_ad[y]>0.) {
      tmp_ad += log(I[y_obj]/Bexp_ad[y]);
      n++;
    }
  }

  if(n>0) q_ad = exp(tmp_ad/n);

  for(y_obj=0,y=ymin_obj-ymin;y_obj<nyr_obj;y_obj++,y++) {
    Ipred_ad[y_obj] = q_ad*Bexp_ad[y];
  }

  //Rprintf("y_obj;\t y;\t I;\t Ipred;\t B;\t Bexp \n");
  //for(y_obj=0,y=ymin_obj-ymin;y_obj<nyr_obj;y_obj++,y++) {
  //  Rprintf("%i;\t %i;\t %f;\t %f;\t %f;\t %f \n",y_obj+ymin_obj,y+ymin,I[y_obj],Ipred_ad[y_obj].getValue(),B_ad[y].getValue(),Bexp_ad[y].getValue());
  //}
  //Rprintf("\n");

}



