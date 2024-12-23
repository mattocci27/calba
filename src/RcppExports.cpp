// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calculate_basal_area_simple
List calculate_basal_area_simple(StringVector sp, NumericVector gx, NumericVector gy, NumericVector ba, double r, bool dist_weighted);
RcppExport SEXP _calba_calculate_basal_area_simple(SEXP spSEXP, SEXP gxSEXP, SEXP gySEXP, SEXP baSEXP, SEXP rSEXP, SEXP dist_weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gx(gxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gy(gySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ba(baSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< bool >::type dist_weighted(dist_weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_basal_area_simple(sp, gx, gy, ba, r, dist_weighted));
    return rcpp_result_gen;
END_RCPP
}
// calculate_basal_area_decay
List calculate_basal_area_decay(NumericVector mu_values, StringVector sp, NumericVector gx, NumericVector gy, NumericVector ba, double r, std::string decay_type);
RcppExport SEXP _calba_calculate_basal_area_decay(SEXP mu_valuesSEXP, SEXP spSEXP, SEXP gxSEXP, SEXP gySEXP, SEXP baSEXP, SEXP rSEXP, SEXP decay_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu_values(mu_valuesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gx(gxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gy(gySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ba(baSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< std::string >::type decay_type(decay_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_basal_area_decay(mu_values, sp, gx, gy, ba, r, decay_type));
    return rcpp_result_gen;
END_RCPP
}
// count_total_cpp
NumericVector count_total_cpp(NumericVector gx, NumericVector gy, double r);
RcppExport SEXP _calba_count_total_cpp(SEXP gxSEXP, SEXP gySEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type gx(gxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gy(gySEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(count_total_cpp(gx, gy, r));
    return rcpp_result_gen;
END_RCPP
}
// count_con_cpp
NumericVector count_con_cpp(StringVector sp, NumericVector gx, NumericVector gy, double r);
RcppExport SEXP _calba_count_con_cpp(SEXP spSEXP, SEXP gxSEXP, SEXP gySEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gx(gxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gy(gySEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(count_con_cpp(sp, gx, gy, r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_calba_calculate_basal_area_simple", (DL_FUNC) &_calba_calculate_basal_area_simple, 6},
    {"_calba_calculate_basal_area_decay", (DL_FUNC) &_calba_calculate_basal_area_decay, 7},
    {"_calba_count_total_cpp", (DL_FUNC) &_calba_count_total_cpp, 3},
    {"_calba_count_con_cpp", (DL_FUNC) &_calba_count_con_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_calba(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
