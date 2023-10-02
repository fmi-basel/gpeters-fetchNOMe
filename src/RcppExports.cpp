// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fetch_data_matrix_from_bams_cpp
Rcpp::List fetch_data_matrix_from_bams_cpp(const Rcpp::CharacterVector& whichContext, const Rcpp::CharacterVector& infiles, const Rcpp::CharacterVector& regionChr, const Rcpp::IntegerVector& regionStart, const Rcpp::IntegerVector& regionEnd, const Rcpp::CharacterVector& seqstring, const Rcpp::IntegerVector& seqStart, const Rcpp::IntegerVector& seqEnd, const Rcpp::LogicalVector& remove_nonunique, const Rcpp::IntegerVector& clip_until_nbg, const Rcpp::NumericVector& max_protect_frac, const Rcpp::NumericVector& max_bisC_meth, const Rcpp::IntegerVector& min_bisC_size, const Rcpp::IntegerVector& mapqMin, const Rcpp::IntegerVector& mapqMax);
RcppExport SEXP _fetchNOMe_fetch_data_matrix_from_bams_cpp(SEXP whichContextSEXP, SEXP infilesSEXP, SEXP regionChrSEXP, SEXP regionStartSEXP, SEXP regionEndSEXP, SEXP seqstringSEXP, SEXP seqStartSEXP, SEXP seqEndSEXP, SEXP remove_nonuniqueSEXP, SEXP clip_until_nbgSEXP, SEXP max_protect_fracSEXP, SEXP max_bisC_methSEXP, SEXP min_bisC_sizeSEXP, SEXP mapqMinSEXP, SEXP mapqMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type whichContext(whichContextSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type infiles(infilesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type regionChr(regionChrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type regionStart(regionStartSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type regionEnd(regionEndSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seqstring(seqstringSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type seqStart(seqStartSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type seqEnd(seqEndSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type remove_nonunique(remove_nonuniqueSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type clip_until_nbg(clip_until_nbgSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type max_protect_frac(max_protect_fracSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type max_bisC_meth(max_bisC_methSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type min_bisC_size(min_bisC_sizeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type mapqMin(mapqMinSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type mapqMax(mapqMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(fetch_data_matrix_from_bams_cpp(whichContext, infiles, regionChr, regionStart, regionEnd, seqstring, seqStart, seqEnd, remove_nonunique, clip_until_nbg, max_protect_frac, max_bisC_meth, min_bisC_size, mapqMin, mapqMax));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fetchNOMe_fetch_data_matrix_from_bams_cpp", (DL_FUNC) &_fetchNOMe_fetch_data_matrix_from_bams_cpp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_fetchNOMe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
