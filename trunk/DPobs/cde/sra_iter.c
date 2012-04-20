
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
void pop_dyn(double,double**,double**,double**);

// global variables
int ymin,ymax,amin,amax,nyr,nag,nit;
double hh,*H,*M;
double *mat,*wght,*sel;
double **C;

extern "C" SEXP run(SEXP R_B0,SEXP R_catch,SEXP R_hh,SEXP R_M,SEXP R_mat,SEXP R_sel,SEXP R_wght,SEXP R_amin,SEXP R_amax,SEXP R_ymin,SEXP R_ymax,SEXP R_nit) {

  //local variables
  int a,y,i,p;

  //////////////
  // Coercion //
  //////////////

  PROTECT(R_B0     = AS_NUMERIC(R_B0));
  PROTECT(R_catch  = AS_NUMERIC(R_catch));
  PROTECT(R_hh     = AS_NUMERIC(R_hh));
  PROTECT(R_M      = AS_NUMERIC(R_M));
  PROTECT(R_mat    = AS_NUMERIC(R_mat));
  PROTECT(R_sel    = AS_NUMERIC(R_sel));
  PROTECT(R_wght   = AS_NUMERIC(R_wght));
  PROTECT(R_amin   = AS_NUMERIC(R_amin));
  PROTECT(R_amax   = AS_NUMERIC(R_amax));
  PROTECT(R_ymin   = AS_NUMERIC(R_ymin));
  PROTECT(R_ymax   = AS_NUMERIC(R_ymax));
  PROTECT(R_nit    = AS_NUMERIC(R_nit));
  
  ////////////////
  // DIMENSIONS //
  ////////////////

  amin = INTEGER_VALUE(R_amin);
  amax = INTEGER_VALUE(R_amax);
  nag  = amax - amin + 1;

  ymin = INTEGER_VALUE(R_ymin);
  ymax = INTEGER_VALUE(R_ymax);
  nyr  = ymax - ymin + 2;
  
  nit  = INTEGER_VALUE(R_nit);
  
  ///////////////////
  // PRELIMINARIES //
  ///////////////////
  
  // B0
  double B0 = NUMERIC_VALUE(R_B0);
 
  // catch
  C = new double*[nyr];
  for(y=0;y<nyr;y++)
    C[y] =  new double[nit];
    
  p = 0;
  for(y=0;y<nyr;y++)
    for(i=0;i<nit;i++) {
      C[y][i] = REAL(R_catch)[p++];
    }

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
  
  double **B    = new double*[nyr];
  double **Bexp = new double*[nyr];
  double **H    = new double*[nyr];
  for(y=0;y<nyr;y++) {
    B[y]    =  new double[nit];
    Bexp[y] =  new double[nit];
    H[y]    =  new double[nit];
  }    
  
  // run model
  pop_dyn(B0,B,Bexp,H);

  ////////////
  // Output //
  ////////////
  int iYear,iIter;
  SEXP dyr,dit,names,dimnames;

  PROTECT(dyr = allocVector(INTSXP,nyr));
  for(iYear=ymin,y=0;y<nyr;iYear++,y++) {
    INTEGER(dyr)[y] = iYear;
  }
  
  PROTECT(dit = allocVector(INTSXP,nit));
  for(iIter=1,i=0;i<nit;iIter++,i++) {
    INTEGER(dit)[i] = iIter;
  }

  // B, Bexp, H
  SEXP _B,_Bexp,_H;
  PROTECT(_B    = allocMatrix(REALSXP,nyr,nit));
  PROTECT(_Bexp = allocMatrix(REALSXP,nyr,nit));
  PROTECT(_H    = allocMatrix(REALSXP,nyr,nit));
  
  p=0;
  for(i=0;i<nit;i++)
    for(y=0;y<nyr;y++) {
      REAL(_B)[p]    = B[y][i];
      REAL(_Bexp)[p] = Bexp[y][i];
      REAL(_H)[p]    = H[y][i];
      p++;
    }

  PROTECT(dimnames = allocVector(VECSXP,2));
  SET_VECTOR_ELT(dimnames,0,dyr);
  SET_VECTOR_ELT(dimnames,1,dit);

  PROTECT(names = allocVector(STRSXP,2));
  SET_STRING_ELT(names,0,mkChar("year"));
  SET_STRING_ELT(names,1,mkChar("iter"));

  setAttrib(dimnames,R_NamesSymbol,names);
  setAttrib(_B,R_DimNamesSymbol,dimnames);
  setAttrib(_Bexp,R_DimNamesSymbol,dimnames);
  setAttrib(_H,R_DimNamesSymbol,dimnames);

  // create combined output list
  SEXP out;
  PROTECT(out = allocVector(VECSXP,3));
  SET_VECTOR_ELT(out,0,_B);
  SET_VECTOR_ELT(out,1,_Bexp);
  SET_VECTOR_ELT(out,2,_H);

  // assign names
  PROTECT(names = allocVector(STRSXP,3));
  SET_STRING_ELT(names,0,mkChar("B"));
  SET_STRING_ELT(names,1,mkChar("Bexp"));
  SET_STRING_ELT(names,2,mkChar("H"));
  setAttrib(out,R_NamesSymbol,names);

  UNPROTECT(22);

  // clean up
  //delete[] B,Bexp,H;

  UNPROTECT(2);
  return out;
  
}

///////////////
// FUNCTIONS //
///////////////

void pop_dyn(double B0,double **B,double **Bexp,double **H) { 

  int a,y,i;
  double alp,bet;

  double rho = 0.;
  double *P = new double[nag];

  double R0;
  double ***N = new double**[nag];
  for(a=0;a<nag;a++) {
    N[a] = new double*[nyr];
    for(y=0;y<nyr;y++) {
      N[a][y] = new double[nit];
    }
  }

  // set up eqm population
  P[0]=1;
  for(a=1;a<nag;a++)
    P[a] = P[a-1]*exp(-M[a-1]);
  P[nag-1] = P[nag-1]/(1-exp(-M[nag-1]));
  for(a=0;a<nag;a++)
    rho += P[a] * mat[a] * wght[a];
  R0 = B0 / rho;
  
  for(i=0;i<nit;i++)
    for(a=0;a<nag;a++)
      N[a][0][i] = R0 * P[a];
  
  // set up S-R parameters

  alp = (4*hh*R0)/(5*hh-1);
  bet = B0*(1-hh)/(5*hh-1);
  
  // Loop through the iters and years
  for(i=0;i<nit;i++) {
  for(y=1;y<nyr;y++) {
  
    // biomass
    B[y-1][i] = 0.;
    Bexp[y-1][i] = 0.;
    for(a=0;a<nag;a++) {
      B[y-1][i] += N[a][y-1][i] * mat[a] * wght[a];
      Bexp[y-1][i] += N[a][y-1][i] * sel[a] * wght[a] * exp(-M[a]/2.);
    }

    // harvest
    H[y-1][i] = 1.;
    if(Bexp[y-1][i]>0.) H[y-1][i]    = C[y-1][i] / Bexp[y-1][i];
    if(H[y-1][i]<0.)    H[y-1][i]    = 0.;
    if(H[y-1][i]>1.)    H[y-1][i]    = 1.;
    if(H[y-1][i]>0.)    Bexp[y-1][i] = C[y-1][i] / H[y-1][i];

    // recruitment

    N[0][y][i] = alp * B[y-1][i]/(bet + B[y-1][i]);

    // adult dynamics

    for(a=1;a<nag;a++)
      N[a][y][i] = N[a-1][y-1][i]*exp(-M[a-1])*(1-sel[a-1]*H[y-1][i]);
    N[nag-1][y][i] = N[nag-1][y][i] + N[nag-1][y-1][i]*exp(-M[nag-1])*(1-sel[nag-1]*H[y-1][i]);

  }
  
  // current biomass
  B[nyr-1][i] = 0.;
  Bexp[nyr-1][i] = 0.;
  for(a=0;a<nag;a++) {
    B[nyr-1][i] += N[a][nyr-1][i] * mat[a] * wght[a];
    Bexp[nyr-1][i] += N[a][nyr-1][i] * sel[a] * wght[a] * exp(-M[a]/2.);
  }
  
  // current H
  H[nyr-1][i] = NA_REAL;
  
  }
  
  // clean up
  //delete[] P;

  //for(a=0;a<nag;a++) {
  //  delete[] N[a];
  //}
  //delete[] N;

}



