// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// concordanceIndex_modified
double concordanceIndex_modified(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, int outx);
RcppExport SEXP _CI_concordanceIndex_modified(SEXP xSEXP, SEXP ySEXP, SEXP deltaXSEXP, SEXP deltaYSEXP, SEXP alphaSEXP, SEXP outxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type deltaX(deltaXSEXP);
    Rcpp::traits::input_parameter< double >::type deltaY(deltaYSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type outx(outxSEXP);
    rcpp_result_gen = Rcpp::wrap(concordanceIndex_modified(x, y, deltaX, deltaY, alpha, outx));
    return rcpp_result_gen;
END_RCPP
}
// shuffle
std::vector<double> shuffle(std::vector<double> array);
RcppExport SEXP _CI_shuffle(SEXP arraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type array(arraySEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle(array));
    return rcpp_result_gen;
END_RCPP
}
// permute_concordanceIndex_modified
std::vector<double> permute_concordanceIndex_modified(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, int outx, int permutations, int nThreads);
RcppExport SEXP _CI_permute_concordanceIndex_modified(SEXP xSEXP, SEXP ySEXP, SEXP deltaXSEXP, SEXP deltaYSEXP, SEXP alphaSEXP, SEXP outxSEXP, SEXP permutationsSEXP, SEXP nThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type deltaX(deltaXSEXP);
    Rcpp::traits::input_parameter< double >::type deltaY(deltaYSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type outx(outxSEXP);
    Rcpp::traits::input_parameter< int >::type permutations(permutationsSEXP);
    Rcpp::traits::input_parameter< int >::type nThreads(nThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(permute_concordanceIndex_modified(x, y, deltaX, deltaY, alpha, outx, permutations, nThreads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _CI_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CI_concordanceIndex_modified", (DL_FUNC) &_CI_concordanceIndex_modified, 6},
    {"_CI_shuffle", (DL_FUNC) &_CI_shuffle, 1},
    {"_CI_permute_concordanceIndex_modified", (DL_FUNC) &_CI_permute_concordanceIndex_modified, 8},
    {"_CI_rcpp_hello", (DL_FUNC) &_CI_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_CI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
