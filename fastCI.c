// #include <Rcpp.h>
// using namespace Rcpp;

// #include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>



void merge_two_sides(double *left_obs, double *left_preds, double *left_disc, double *left_pairs, int left_length,
                     double *right_obs, double *right_preds, double *right_disc, double *right_pairs, int right_length,
                     double *out_obs, double *out_preds, double *out_disc, double *out_pairs, int out_length,
                     int discard_x_ties, int discard_y_ties){

  //RLR = Right List Removed  
  int RLR = 0;
  //LLL = Left List Left
  int LLL = left_length;
  
  // Length Right
  int LR = right_length;

  //Left Index; Right Index, index (of output vector)
  int Li = 0;
  int Ri = 0;
  int i = 0;

  // Handle any potential ties in the predictions across L and R, which as R_pred >= L_pred, 
  // can only emerge from having elements in the right side equal to max(L_pred)
  

  while(i < out_length){
    
    if(LLL == 0){
      //// If left list is empty the only things we can do is fill in the
      //// output with right list elements.
      out_obs[i] = right_obs[Ri];
      out_preds[i] = right_preds[Ri];
      out_disc[i] = right_disc[Ri] + LLL; //LLL = 0, but for consistency leaving here;
      out_pairs[i] = right_pairs[Ri];
      Ri = Ri + 1;
      i = i + 1;
      continue;
    }
    if(RLR == LR){
      //// If all elements from the right list have been removed, we fill in from left list
      out_obs[i] = left_obs[Li];
      out_preds[i] = left_preds[Li];
      out_disc[i] = left_disc[Li] + RLR;
      out_pairs[i] = left_pairs[Li];
      Li = Li + 1;
      i = i + 1;
      continue;
    }
    // if(left_preds[Li] == right_preds[Ri] || (left_obs[Li] == right_obs[Ri] && outx)){
    if(left_preds[Li] == right_preds[Ri] || (left_obs[Li] == right_obs[Ri])){

      // Is this still wrong?
      //// This loop removes elements from the left list while they remain tied with the leftmost element of the right list 
      while(LLL && (left_obs[Li] == right_obs[Ri] || left_preds[Li] == right_preds[Ri])){
        out_obs[i] = left_obs[Li];
        out_preds[i] = left_preds[Li];
        out_disc[i] = left_disc[Li] + RLR;
        out_pairs[i] = left_pairs[Li] - 1;
        i = i + 1;
        LLL = LLL - 1;
        Li = Li + 1;
      }
      out_obs[i] = right_obs[Ri];
      out_preds[i] = right_preds[Ri];
      out_disc[i] = right_disc[Ri] + LLL ;
      out_pairs[i] = right_pairs[Ri] - 1;
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    } else if(left_obs[Li] < right_obs[Ri]) {
      out_obs[i] = left_obs[Li];
      out_preds[i] = left_preds[Li];
      out_disc[i] = left_disc[Li] + RLR;
      out_pairs[i] = left_pairs[Li];
      LLL = LLL - 1;
      Li = Li + 1;
      i = i + 1;
    } else if(left_obs[Li] > right_obs[Ri]) {
      out_obs[i] = right_obs[Ri];
      out_preds[i] = right_preds[Ri];
      out_disc[i] = right_disc[Ri] + LLL;
      out_pairs[i] = right_pairs[Ri];
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    } else {
      // only case left is if the two values are equal and outx is false
      // stop("Not implemented correctly?")
      out_obs[i] = left_obs[Li];
      out_preds[i] = left_preds[Li];
      out_disc[i] = left_disc[Li] + RLR + 0.5;
      out_pairs[i] = left_pairs[Li];
      i = i + 1;
      out_obs[i] = right_obs[Ri];
      out_preds[i] = right_preds[Ri];
      out_disc[i] = right_disc[Ri] + LLL - 0.5;
      out_pairs[i] = right_pairs[Ri];
      LLL = LLL - 1;
      Li = Li + 1;
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    }
  }
  return;
}


void merge_sort(double *in_obs, double *in_preds, double *in_disc, double *in_pairs, int in_length, 
                double *out_obs, double *out_preds, double *out_disc, double *out_pairs, int out_length,
                int discard_x_ties, int discard_y_ties){
  if(in_length == 1){
    return;
  } else {

    int split_idx = floor(in_length/2);

    double *left_obs = &in_obs[0];
    double *left_preds = &in_preds[0];
    double *left_disc = &in_disc[0];
    double *left_pairs = &in_pairs[0];

    double *oleft_obs = &out_obs[0];
    double *oleft_preds = &out_preds[0];
    double *oleft_disc = &out_disc[0];
    double *oleft_pairs = &out_pairs[0];

    int left_length = split_idx;

    double *right_obs = &in_obs[split_idx];
    double *right_preds = &in_preds[split_idx];
    double *right_disc = &in_disc[split_idx];
    double *right_pairs = &in_pairs[split_idx];
    
    double *oright_obs = &out_obs[split_idx];
    double *oright_preds = &out_preds[split_idx];
    double *oright_disc = &out_disc[split_idx];
    double *oright_pairs = &out_pairs[split_idx];

    int right_length = in_length - split_idx;
    

    merge_sort(oleft_obs, oleft_preds, oleft_disc, oleft_pairs, left_length, 
               left_obs, left_preds, left_disc, left_pairs, left_length,
               discard_x_ties, discard_y_ties);
    merge_sort(oright_obs, oright_preds, oright_disc, oright_pairs, right_length, 
               right_obs, right_preds, right_disc, right_pairs, right_length,
               discard_x_ties, discard_y_ties);
    merge_two_sides(left_obs, left_preds, left_disc, left_pairs, left_length,
                    right_obs, right_preds, right_disc, right_pairs, right_length,
                    out_obs, out_preds, out_disc, out_pairs, out_length,
                    discard_x_ties, discard_y_ties);
    return;
  }
}

SEXP merge_sort_c(SEXP pin_obs,
                  SEXP pin_preds,
                  SEXP pin_disc,
                  SEXP pin_pairs,
                  SEXP pn,
                  SEXP pdiscard_x_ties, 
                  SEXP pdiscard_y_ties){


  int N = *INTEGER(pn);
  // Problem here
  int discard_x_ties = *INTEGER(pdiscard_x_ties);
  int discard_y_ties = *INTEGER(pdiscard_y_ties);


  double *in_obs = REAL(PROTECT(allocVector(REALSXP, N)));
  double *in_preds = REAL(PROTECT(allocVector(REALSXP, N)));
  double *in_disc = REAL(PROTECT(allocVector(REALSXP, N)));
  double *in_pairs = REAL(PROTECT(allocVector(REALSXP, N)));

  SEXP pout_obs = PROTECT(allocVector(REALSXP, N));
  SEXP pout_preds = PROTECT(allocVector(REALSXP, N));
  SEXP pout_disc = PROTECT(allocVector(REALSXP, N));
  SEXP pout_pairs = PROTECT(allocVector(REALSXP, N));
  SEXP out = PROTECT(allocVector(VECSXP, 4));

  double *out_obs = REAL(pout_obs);
  double *out_preds = REAL(pout_preds);
  double *out_disc = REAL(pout_disc);
  double *out_pairs = REAL(pout_pairs);

  memcpy(out_obs, REAL(pin_obs), N * sizeof(*out_obs));
  memcpy(out_preds, REAL(pin_preds), N * sizeof(*out_preds));
  memcpy(out_disc, REAL(pin_disc), N * sizeof(*out_disc));
  memcpy(out_pairs, REAL(pin_pairs), N * sizeof(*out_pairs));
  
  memcpy(in_obs, REAL(pin_obs), N * sizeof(*in_obs));
  memcpy(in_preds, REAL(pin_preds), N * sizeof(*in_preds));
  memcpy(in_disc, REAL(pin_disc), N * sizeof(*in_disc));
  memcpy(in_pairs, REAL(pin_pairs), N * sizeof(*in_pairs));


  merge_sort(in_obs, in_preds, in_disc, in_pairs, N, 
             out_obs, out_preds, out_disc, out_pairs, N,
             discard_x_ties, discard_y_ties);

  SET_VECTOR_ELT(out, 0, pout_obs);
  SET_VECTOR_ELT(out, 1, pout_preds);
  SET_VECTOR_ELT(out, 2, pout_disc);
  SET_VECTOR_ELT(out, 3, pout_pairs);

  UNPROTECT(9);
  return out;

}


static const R_CallMethodDef callMethods[]  = {
  {"merge_sort_c", (DL_FUNC) &merge_sort_c, 6},
  {NULL, NULL, 0}
};

// 
// int main(){

//   double ci = 0;
//   int outx = 1;
//   int n = 6;

//   double obs[] = {1,2,3,4,6,5};
//   double preds[] = {1,2,3,4,5,5};
//   double disc[] = {0,0,0,0,0,0};
//   double pairs[] = {5,5,5,5,5,5};

//   double obso[n];
//   double predso[n];
//   double disco[n];
//   double pairso[n];

//   for (int i = 0; i < n; i++){
//     obso[i] = obs[i];
//     predso[i] = preds[i];
//     disco[i] = disc[i];
//     pairso[i] = pairs[i];
//   }

//   merge_sort(obs, preds, disc, pairs, n, 
//                obso, predso, disco, pairso, n,
//                outx);

//   for (int i = 0; i < n; i++){
//     ci += disco[i];
//   }
//   ci = ci/2.0;

//   printf("Inverted Pairs = %f", ci);

// }







