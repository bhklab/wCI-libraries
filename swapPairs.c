// #include <Rcpp.h>
// using namespace Rcpp;

// #include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

#include "xoroshiro128+.h"



void swap_consecutive(double *vector, int index){

	double temp = vector[index];
	vector[index] = vector[index+1];
	vector[index+1] = temp;

} 

SEXP do_swaps_c(SEXP pvector,
                  SEXP pvlen,
                  SEXP pnswap,
                  SEXP pseed){


  double vector_length = *REAL(pvlen);

  double nswap = *REAL(pnswap);
  SEXP svector = PROTECT(allocVector(REALSXP, vector_length));

  double *vector = REAL(svector);
  double *seed = REAL(pseed);


  memcpy(vector, REAL(pvector), (uint64_t)vector_length * sizeof(*vector));
  // char *statechar = (char []) {0x77, 0x87, 0xfe, 0xb6, 0x8e, 0xe0, 0x03, 0x73, 0xab, 0x2d, 0xfa, 0x72, 0x8c, 0xa0, 0x5c, 0x26};
  uint64_t *state = (uint64_t*) seed;

  uint64_t index;

  for(double i = 0; i < nswap; i++){
  	index = generate_random_index(state, vector_length - 1);
  	swap_consecutive(vector, index);
  }


  UNPROTECT(1);
  return svector;

}


static const R_CallMethodDef callMethods[]  = {
  {"do_swaps_c", (DL_FUNC) &do_swaps_c, 3},
  {NULL, NULL, 0}
};

// int main(){

// 	char *statechar = (char []) {0x77, 0x87, 0xfe, 0xb6, 0x8e, 0xe0, 0x03, 0x73, 0xab, 0x2d, 0xfa, 0x72, 0x8c, 0xa0, 0x5c, 0x26};
//   uint64_t *state = (uint64_t*) statechar;

//   uint64_t number;
//   uint64_t maxindex = 10000;

//   number = generate_random_index(state, maxindex);

// 	printf("Random Number = %llu \n", number);
//   number = generate_random_index(state, maxindex);
// 	printf("Random Number = %llu \n", number);


// }
