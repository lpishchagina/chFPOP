// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_changepoints
List get_changepoints(Rcpp::NumericMatrix data, double penalty, std::string method, bool showNbCands);
RcppExport SEXP _chFPOP_get_changepoints(SEXP dataSEXP, SEXP penaltySEXP, SEXP methodSEXP, SEXP showNbCandsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type showNbCands(showNbCandsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_changepoints(data, penalty, method, showNbCands));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chFPOP_get_changepoints", (DL_FUNC) &_chFPOP_get_changepoints, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_chFPOP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
