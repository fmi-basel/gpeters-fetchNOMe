// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fetch_molec_data_list_from_bams_cpp
Rcpp::List fetch_molec_data_list_from_bams_cpp(const Rcpp::CharacterVector& whichContext, const Rcpp::CharacterVector& infiles, const Rcpp::CharacterVector& regionChr, const Rcpp::IntegerVector& regionStart, const Rcpp::IntegerVector& regionEnd, const Rcpp::CharacterVector& seqstring, const Rcpp::IntegerVector& seqStart, const Rcpp::IntegerVector& seqEnd, const Rcpp::LogicalVector& remove_nonunique, const Rcpp::IntegerVector& clip_until_nbg, const Rcpp::NumericVector& max_protect_frac, const Rcpp::NumericVector& max_bisC_meth, const Rcpp::IntegerVector& min_bisC_size, const Rcpp::IntegerVector& mapqMin, const Rcpp::IntegerVector& mapqMax, const Rcpp::CharacterVector& alignerUsed, const Rcpp::CharacterVector& SMFenzymeUsed, const Rcpp::IntegerVector& min_frag_data_len, const Rcpp::NumericVector& min_frag_data_dens, const Rcpp::LogicalVector& data_as_rle);
RcppExport SEXP _fetchNOMe_fetch_molec_data_list_from_bams_cpp(SEXP whichContextSEXP, SEXP infilesSEXP, SEXP regionChrSEXP, SEXP regionStartSEXP, SEXP regionEndSEXP, SEXP seqstringSEXP, SEXP seqStartSEXP, SEXP seqEndSEXP, SEXP remove_nonuniqueSEXP, SEXP clip_until_nbgSEXP, SEXP max_protect_fracSEXP, SEXP max_bisC_methSEXP, SEXP min_bisC_sizeSEXP, SEXP mapqMinSEXP, SEXP mapqMaxSEXP, SEXP alignerUsedSEXP, SEXP SMFenzymeUsedSEXP, SEXP min_frag_data_lenSEXP, SEXP min_frag_data_densSEXP, SEXP data_as_rleSEXP) {
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
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type alignerUsed(alignerUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type SMFenzymeUsed(SMFenzymeUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type min_frag_data_len(min_frag_data_lenSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type min_frag_data_dens(min_frag_data_densSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type data_as_rle(data_as_rleSEXP);
    rcpp_result_gen = Rcpp::wrap(fetch_molec_data_list_from_bams_cpp(whichContext, infiles, regionChr, regionStart, regionEnd, seqstring, seqStart, seqEnd, remove_nonunique, clip_until_nbg, max_protect_frac, max_bisC_meth, min_bisC_size, mapqMin, mapqMax, alignerUsed, SMFenzymeUsed, min_frag_data_len, min_frag_data_dens, data_as_rle));
    return rcpp_result_gen;
END_RCPP
}
// fetch_cooc_ctable_from_bams_cpp
Rcpp::List fetch_cooc_ctable_from_bams_cpp(const Rcpp::CharacterVector& infiles, const Rcpp::CharacterVector& regionChr, const Rcpp::IntegerVector& regionStart, const Rcpp::IntegerVector& regionEnd, const Rcpp::IntegerVector& max_spac, const Rcpp::CharacterVector& seqstring, const Rcpp::IntegerVector& seqStart, const Rcpp::IntegerVector& seqEnd, const Rcpp::LogicalVector& remove_nonunique, const Rcpp::IntegerVector& clip_until_nbg, const Rcpp::NumericVector& max_protect_frac, const Rcpp::NumericVector& max_bisC_meth, const Rcpp::IntegerVector& min_bisC_size, const Rcpp::IntegerVector& mapqMin, const Rcpp::IntegerVector& mapqMax, const Rcpp::CharacterVector& alignerUsed, const Rcpp::CharacterVector& SMFenzymeUsed, const Rcpp::IntegerVector& min_frag_data_len, const Rcpp::NumericVector& min_frag_data_dens);
RcppExport SEXP _fetchNOMe_fetch_cooc_ctable_from_bams_cpp(SEXP infilesSEXP, SEXP regionChrSEXP, SEXP regionStartSEXP, SEXP regionEndSEXP, SEXP max_spacSEXP, SEXP seqstringSEXP, SEXP seqStartSEXP, SEXP seqEndSEXP, SEXP remove_nonuniqueSEXP, SEXP clip_until_nbgSEXP, SEXP max_protect_fracSEXP, SEXP max_bisC_methSEXP, SEXP min_bisC_sizeSEXP, SEXP mapqMinSEXP, SEXP mapqMaxSEXP, SEXP alignerUsedSEXP, SEXP SMFenzymeUsedSEXP, SEXP min_frag_data_lenSEXP, SEXP min_frag_data_densSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type infiles(infilesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type regionChr(regionChrSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type regionStart(regionStartSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type regionEnd(regionEndSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type max_spac(max_spacSEXP);
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
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type alignerUsed(alignerUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type SMFenzymeUsed(SMFenzymeUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type min_frag_data_len(min_frag_data_lenSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type min_frag_data_dens(min_frag_data_densSEXP);
    rcpp_result_gen = Rcpp::wrap(fetch_cooc_ctable_from_bams_cpp(infiles, regionChr, regionStart, regionEnd, max_spac, seqstring, seqStart, seqEnd, remove_nonunique, clip_until_nbg, max_protect_frac, max_bisC_meth, min_bisC_size, mapqMin, mapqMax, alignerUsed, SMFenzymeUsed, min_frag_data_len, min_frag_data_dens));
    return rcpp_result_gen;
END_RCPP
}
// fetch_data_matrix_from_bams_cpp
Rcpp::List fetch_data_matrix_from_bams_cpp(const Rcpp::CharacterVector& whichContext, const Rcpp::CharacterVector& infiles, const Rcpp::CharacterVector& regionChr, const Rcpp::IntegerVector& regionStart, const Rcpp::IntegerVector& regionEnd, const Rcpp::CharacterVector& seqstring, const Rcpp::IntegerVector& seqStart, const Rcpp::IntegerVector& seqEnd, const Rcpp::LogicalVector& remove_nonunique, const Rcpp::IntegerVector& clip_until_nbg, const Rcpp::NumericVector& max_protect_frac, const Rcpp::NumericVector& max_bisC_meth, const Rcpp::IntegerVector& min_bisC_size, const Rcpp::IntegerVector& mapqMin, const Rcpp::IntegerVector& mapqMax, const Rcpp::CharacterVector& alignerUsed, const Rcpp::CharacterVector& SMFenzymeUsed, const Rcpp::IntegerVector& min_frag_data_len, const Rcpp::NumericVector& min_frag_data_dens);
RcppExport SEXP _fetchNOMe_fetch_data_matrix_from_bams_cpp(SEXP whichContextSEXP, SEXP infilesSEXP, SEXP regionChrSEXP, SEXP regionStartSEXP, SEXP regionEndSEXP, SEXP seqstringSEXP, SEXP seqStartSEXP, SEXP seqEndSEXP, SEXP remove_nonuniqueSEXP, SEXP clip_until_nbgSEXP, SEXP max_protect_fracSEXP, SEXP max_bisC_methSEXP, SEXP min_bisC_sizeSEXP, SEXP mapqMinSEXP, SEXP mapqMaxSEXP, SEXP alignerUsedSEXP, SEXP SMFenzymeUsedSEXP, SEXP min_frag_data_lenSEXP, SEXP min_frag_data_densSEXP) {
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
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type alignerUsed(alignerUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type SMFenzymeUsed(SMFenzymeUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type min_frag_data_len(min_frag_data_lenSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type min_frag_data_dens(min_frag_data_densSEXP);
    rcpp_result_gen = Rcpp::wrap(fetch_data_matrix_from_bams_cpp(whichContext, infiles, regionChr, regionStart, regionEnd, seqstring, seqStart, seqEnd, remove_nonunique, clip_until_nbg, max_protect_frac, max_bisC_meth, min_bisC_size, mapqMin, mapqMax, alignerUsed, SMFenzymeUsed, min_frag_data_len, min_frag_data_dens));
    return rcpp_result_gen;
END_RCPP
}
// fetch_protect_stats_from_bams_cpp
Rcpp::List fetch_protect_stats_from_bams_cpp(const Rcpp::CharacterVector& infiles, const Rcpp::CharacterVector& regionChr, const Rcpp::IntegerVector& regionStart, const Rcpp::IntegerVector& regionEnd, const Rcpp::CharacterVector& seqstring, const Rcpp::IntegerVector& seqStart, const Rcpp::IntegerVector& seqEnd, const Rcpp::LogicalVector& remove_nonunique, const Rcpp::IntegerVector& clip_until_nbg, const Rcpp::NumericVector& max_protect_frac, const Rcpp::NumericVector& max_bisC_meth, const Rcpp::IntegerVector& min_bisC_size, const Rcpp::IntegerVector& mapqMin, const Rcpp::IntegerVector& mapqMax, const Rcpp::CharacterVector& alignerUsed, const Rcpp::CharacterVector& SMFenzymeUsed, const Rcpp::IntegerVector& min_frag_data_len, const Rcpp::NumericVector& min_frag_data_dens);
RcppExport SEXP _fetchNOMe_fetch_protect_stats_from_bams_cpp(SEXP infilesSEXP, SEXP regionChrSEXP, SEXP regionStartSEXP, SEXP regionEndSEXP, SEXP seqstringSEXP, SEXP seqStartSEXP, SEXP seqEndSEXP, SEXP remove_nonuniqueSEXP, SEXP clip_until_nbgSEXP, SEXP max_protect_fracSEXP, SEXP max_bisC_methSEXP, SEXP min_bisC_sizeSEXP, SEXP mapqMinSEXP, SEXP mapqMaxSEXP, SEXP alignerUsedSEXP, SEXP SMFenzymeUsedSEXP, SEXP min_frag_data_lenSEXP, SEXP min_frag_data_densSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
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
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type alignerUsed(alignerUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type SMFenzymeUsed(SMFenzymeUsedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type min_frag_data_len(min_frag_data_lenSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type min_frag_data_dens(min_frag_data_densSEXP);
    rcpp_result_gen = Rcpp::wrap(fetch_protect_stats_from_bams_cpp(infiles, regionChr, regionStart, regionEnd, seqstring, seqStart, seqEnd, remove_nonunique, clip_until_nbg, max_protect_frac, max_bisC_meth, min_bisC_size, mapqMin, mapqMax, alignerUsed, SMFenzymeUsed, min_frag_data_len, min_frag_data_dens));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fetchNOMe_fetch_molec_data_list_from_bams_cpp", (DL_FUNC) &_fetchNOMe_fetch_molec_data_list_from_bams_cpp, 20},
    {"_fetchNOMe_fetch_cooc_ctable_from_bams_cpp", (DL_FUNC) &_fetchNOMe_fetch_cooc_ctable_from_bams_cpp, 19},
    {"_fetchNOMe_fetch_data_matrix_from_bams_cpp", (DL_FUNC) &_fetchNOMe_fetch_data_matrix_from_bams_cpp, 19},
    {"_fetchNOMe_fetch_protect_stats_from_bams_cpp", (DL_FUNC) &_fetchNOMe_fetch_protect_stats_from_bams_cpp, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_fetchNOMe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
