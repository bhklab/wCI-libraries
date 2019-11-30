/*
// Fast permutations for rCI using a naive matrix based approach.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

#include "xoroshiro128+.h"


struct runningStat{
    double numerator;
    double denominator;
};

typedef struct runningStat runningStat;

void printVec(uint64_t *list, uint64_t N){
  for(uint64_t i = 0; i < N; i ++){
    printf("%lld, ", list[i]);
  }
}

// Using "inside out" Fisher Yates https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
void sampleIdx(uint64_t N, uint64_t *permPointer, uint64_t *state){
  uint64_t j;
  permPointer[0] = 0;

  for(uint64_t i = 1; i <= N-1; i++){
      j = generate_random_index(state, i);
      if(j != i){
        permPointer[i] = permPointer[j];
      }
      permPointer[j] = i;
  }

};
  
  


runningStat returnCurrentIncrement(uint64_t xVal, uint64_t yVal, int xties, int yties){

  runningStat res;

  if(xties == yties & yties == 0){
    // printf("Ignore ties\n. ");
      res.numerator = (double) ((xVal == yVal) & (xVal != 0) & (yVal != 0)); // Datasets agree and are non-zero
      res.denominator = (double) ((xVal != 0) & (yVal != 0)); // Both are non-zero
    } else if (xties == yties & yties == 1){

      res.numerator = ((double) xVal == yVal)*0.5 + (double) ((xVal != 0) & (yVal != 0)) * 0.5 ; // Datasets agree (0.5) non zero (0.5)
      res.denominator = 1; // All pairs are counted. 

    } else {
      if (xties == 1 & yties == 0){
        res.numerator = (((double) xVal == yVal)*0.5 + (double) (xVal != 0) * 0.5) * (double) (yVal != 0);
        res.denominator = (double) (yVal != 0);
      }
      if (xties == 0 & yties == 1){
        res.numerator = (((double) xVal == yVal)*0.5 + (double) (yVal != 0) * 0.5) * (double) (xVal != 0);
        res.denominator = (double) (xVal != 0);
      }
    }
    return res;
};


uint64_t returnLinearIdx(uint64_t row, uint64_t column, uint64_t N){
  return(column*N + row);
};



double rciBoot(int *xmat, int *ymat, double obsCI, uint64_t R, uint64_t B, uint64_t N, int xties, int yties, uint64_t *state){

  double currCI;
  uint64_t totalSeenLarger = 0;
  uint64_t i = 0;
  // uint64_t j = 0;
  uint64_t *permIdx = malloc(N * sizeof(uint64_t));


  uint64_t xVal;
  uint64_t yVal;

  runningStat RS;
  RS.numerator = 0;
  RS.denominator = 0;
  runningStat CS;

  while(i < B){

    sampleIdx(N, permIdx, state);
    // printVec(permIdx, N);
    // printf("\n");

    for(uint64_t j = 0; j < N; j++){

      for(uint64_t k = 0; k < N; k++){

        xVal = xmat[returnLinearIdx(j,k,N)];
        yVal = ymat[returnLinearIdx(permIdx[j],permIdx[k],N)];
        // printf("%ld,%ld. ", xVal, yVal);

        CS = returnCurrentIncrement(xVal, yVal, xties, yties);
        // printf("CurNum: %f. \n CurDenom: %f.\n", CS.numerator, CS.denominator);

        RS.numerator += CS.numerator;
        RS.denominator += CS.denominator;

      }

    }

    currCI = ((double)RS.numerator)/((double)RS.denominator);
    // printf("Obs Num: %f. \n Obs Denom: %f.\n", RS.numerator, RS.denominator);

    // printf("Obs CI: %f. \n Curr CI: %f.\n", obsCI, currCI);
    if(fabs(currCI - 0.5) >= fabs(obsCI - 0.5)){
      totalSeenLarger++;
    }

    if(totalSeenLarger==R){
      // printf("i is: %ld. Escape from loop\n", i);
      break;
    }
    RS.numerator = 0;
    RS.denominator = 0;
    i++;
  }

  // TODO:: Figure out how the 0 based indexing affects calculation of the pvalues
  // printf("%ld\n", totalSeenLarger);

  free(permIdx);

  if(totalSeenLarger==R){
    return(((double)totalSeenLarger)/((double)i+1));
  }
  return(((double)(totalSeenLarger + 1))/((double)(B + 1)));

}


SEXP permC(SEXP pin_x,
             SEXP pin_y,
             SEXP pobsCI,
             SEXP pR,
             SEXP pB,
             SEXP pn,
             SEXP pxties,
             SEXP pyties, 
             SEXP pseed){
  
  double Ndouble = *REAL(pn);
  
  // int discard_x_ties = *INTEGER(pdiscard_x_ties);
  // int discard_y_ties = *INTEGER(pdiscard_y_ties);
  
  double obsCI = *REAL(pobsCI);

  double Rdouble = *REAL(pR);  
  double Bdouble = *REAL(pB);

  uint64_t N = (uint64_t) Ndouble;
  uint64_t R = (uint64_t) Rdouble;
  uint64_t B = (uint64_t) Bdouble;

  int xties = *INTEGER(pxties);
  int yties = *INTEGER(pyties);
  
  SEXP pout = PROTECT(allocVector(REALSXP,1));
  
  double *out = REAL(pout);
  
  double *seed = REAL(pseed);
  uint64_t *state = (uint64_t*) seed;

  
  out[0] = rciBoot(INTEGER(pin_x), INTEGER(pin_y), obsCI, R, B, N, xties, yties, state);
  // printf("%f\n", out[0]);
  
  UNPROTECT(1);
  
  return pout;
  
}


