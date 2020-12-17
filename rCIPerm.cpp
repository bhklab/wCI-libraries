#include <Rcpp.h>

extern "C" SEXP permC(SEXP pin_x,
             SEXP pin_y,
             SEXP pobsCI,
             SEXP pR,
             SEXP pB,
             SEXP pn,
             SEXP pxties,
             SEXP pyties, 
             SEXP palternative, 
             SEXP pseed);

extern "C" SEXP newPairedConcIndex(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pdeltaX,
             SEXP pdeltaY,
             SEXP plogic);


extern "C" SEXP bootC(SEXP prcimat,
             SEXP pR,
             SEXP pn,
             SEXP pxties,
             SEXP pyties, 
             SEXP pseed);

extern "C" SEXP kernalizedPairedConcIndex(SEXP pin_x,
                                SEXP pin_y,
                                SEXP pn,
                                SEXP psymmetric,
                                SEXP pkern1,
                                SEXP pkern2);

extern "C" SEXP rCI(SEXP pin_x,
                  SEXP pin_y,
                  SEXP pdelta_x,
                  SEXP pdelta_y,
                  SEXP pn,
                  SEXP pdiscard_x_ties, 
                  SEXP pdiscard_y_ties);


// [[Rcpp::export]]
SEXP rCIPermC(SEXP pin_x, SEXP pin_y, SEXP pobsCI, SEXP pR, SEXP pB, SEXP pn, SEXP pxties, SEXP pyties, SEXP palternative,  SEXP pseed) {
    return permC(pin_x,
             pin_y,
             pobsCI,
             pR,
             pB,
             pn,
             pxties,
             pyties,
             palternative, 
             pseed);
}


// [[Rcpp::export]]
SEXP rCIBootC(SEXP prcimat, SEXP pR, SEXP pn, SEXP pxties, SEXP pyties, SEXP pseed) {
    return bootC(prcimat,
             pR,
             pn,
             pxties,
             pyties,
             pseed);
}

// [[Rcpp::export]]
SEXP newPCI(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pdeltaX,
             SEXP pdeltaY,
             SEXP plogic) {
    return newPairedConcIndex(pin_x,
             pin_y,
             pn,
             pxties,
             pyties,
             pdeltaX,
             pdeltaY,
             plogic);
}


// [[Rcpp::export]]
SEXP frCI(SEXP pin_x,
                  SEXP pin_y,
                  SEXP pdelta_x,
                  SEXP pdelta_y,
                  SEXP pn,
                  SEXP pdiscard_x_ties, 
                  SEXP pdiscard_y_ties){

        return rCI(pin_x,
                  pin_y,
                  pdelta_x,
                  pdelta_y,
                  pn,
                  pdiscard_x_ties, 
                  pdiscard_y_ties);
}

// [[Rcpp::export]]
SEXP KCI(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP psymmetric,
             SEXP pkern1,
             SEXP pkern2) {
    return kernalizedPairedConcIndex(pin_x,
             pin_y,
             pn,
             psymmetric,
             pkern1,
             pkern2);
}
