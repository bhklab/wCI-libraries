// SEXP rCIBoot(SEXP pin_x,
//          SEXP pin_y,
//          SEXP pdelta_x,
//          SEXP pdelta_y,
//          SEXP pn,
//          SEXP pdiscard_x_ties, 
//          SEXP pdiscard_y_ties,
//          SEXP pR){
  
//     long N = *INTEGER(pn);
  
//     // int discard_x_ties = *INTEGER(pdiscard_x_ties);
//     // int discard_y_ties = *INTEGER(pdiscard_y_ties);
    
//     double delta_x = *REAL(pdelta_x);
//     double delta_y = *REAL(pdelta_y);
    
//     double *in_x = REAL(PROTECT(allocVector(REALSXP, N)));
//     double *in_y = REAL(PROTECT(allocVector(REALSXP, N)));
    
//     long R = *INTEGER(pR);
    
//     SEXP pout_t0 = PROTECT(allocVector(REALSXP,1));
//     SEXP pout_t = PROTECT(allocVector(REALSXP,R));
    
//     double *out_rci = REAL(pout_rci);
//     double *out_valid = REAL(pout_valid);
    
//     memcpy(in_x, REAL(pin_x), N * sizeof(*in_x));
//     memcpy(in_y, REAL(pin_y), N * sizeof(*in_y));
    
    
    
//     SEXP out = PROTECT(allocVector(VECSXP, 2));
    
//     // printf("First x element is: %f\n", in_x[0]);
//     Res res = rci(in_x, in_y, delta_x, delta_y, N);
//     out_rci[0] = res.rCI;
//     out_valid[0] = res.valid_pairs;
    
//     // TODO: these should be REALSXP pointers
//     SET_VECTOR_ELT(out, 0, pout_t0);
//     SET_VECTOR_ELT(out, 1, pout_t);
    
//     UNPROTECT(5);
    
//     return out;
    
// }

// static const R_CallMethodDef callMethods[]  = {
//   {"rCIBoot", (DL_FUNC) &rCIBoot, 7},
//   {NULL, NULL, 0}
// };

