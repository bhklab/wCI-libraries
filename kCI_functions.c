#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

#include "xoroshiro128+.h"



struct retRes{
    double C;
    double D;
    double CC;
    double CD;
    double DD;
    uint64_t N;
};

typedef struct retRes retRes;



double bKernel(double diff, double kern1, double kern2){

  return(1.0/(1.0+exp(kern2*(fabs(diff) - kern2))));

}

void computePairKCI(double x1, double x2, double y1, double y2, 
  double *C1, double *C2, double *D1, double *D2, int symmetric, double kern1, double kern2){

  double diffX = x1 - x2;
  double diffY = y1 - y2;
  double weight;

  if(symmetric == 0){
    weight = bKernel(diffX, kern1, kern2);
  } else {
    weight = bKernel(diffX, kern1, kern2)*bKernel(diffY, kern1, kern2);

  }

  if((diffX > 0 & diffY > 0) | (diffX < 0 & diffY < 0)){

    *C1 += weight;
    *C2 += weight;    

  } else if((diffX > 0 & diffY < 0) | (diffX < 0 & diffY > 0)) {

    *D1 += weight;
    *D2 += weight;
  
  }

  return;
}


retRes kernalized_concordance_index(double *x, double *y, uint64_t N, int symmetric, double kern1, double kern2) {

  double c[N];
  double d[N];

  // std::list<bool> cdseq;

  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }


  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {

      computePairKCI(x[i], x[j], y[i], y[j], &c[i], &c[j], &d[i], &d[j], symmetric, kern1, kern2);

    }
  }

  double C = 0.0;
  double D = 0.0;
  double CC = 0.0;
  double CD = 0.0;
  double DD = 0.0;



  for (int i = 0; i < N; ++i) {
    C += c[i];
    D += d[i];
    CC += c[i] * (c[i] - 1);
    DD += d[i] * (d[i] - 1);
    CD += c[i] * d[i];
  }


  retRes ret;
  ret.C = C;
  ret.D = D;
  ret.CC = CC;
  ret.DD = DD;
  ret.CD = CD;
  ret.N = N;
  // reteq.cdseq = cdseq;
  return ret;


}



SEXP kernalizedPairedConcIndex(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP psymmetric,
             SEXP pkern1,
             SEXP pkern2){
  
  double Ndouble = *REAL(pn);
  
  // int discard_x_ties = *INTEGER(pdiscard_x_ties);
  // int discard_y_ties = *INTEGER(pdiscard_y_ties);

  uint64_t N = (uint64_t) Ndouble;
  
  
  SEXP pout = PROTECT(allocVector(REALSXP,6));
  
  double *out = REAL(pout);
  
  
  retRes res = kernalized_concordance_index(REAL(pin_x), 
                                            REAL(pin_y), 
                                            N, 
                                            *INTEGER(psymmetric),
                                            *REAL(pkern1),
                                            *REAL(pkern2));
  

  out[0] = res.C;
  out[1] = res.D; 
  out[2] = res.CC;
  out[3] = res.DD;
  out[4] = res.CD;
  out[5] = res.N;

  // printf("%f\n", out[0]);
  
  UNPROTECT(1);
  
  return pout;
  
}
