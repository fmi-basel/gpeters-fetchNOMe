#ifndef _functions_export_rcpp_hpp_
#define _functions_export_rcpp_hpp_


#include <Rcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <vector>
#include <map>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdint>
#include <stdbool.h>
#include "refSeqInfo-class.h"
#include "regionData-class.h"
#include "coocCntTable-class.h"
#include "protectStats-class.h"
#include "fetch_data_from_bam.h"


using namespace std;


bool _VERBOSE_ = 0;

// [[Rcpp::export]]
Rcpp::List fetch_molec_data_list_from_bams_cpp(const Rcpp::CharacterVector& whichContext,
                                               const Rcpp::CharacterVector& infiles,
                                               const Rcpp::CharacterVector& regionChr,
                                               const Rcpp::IntegerVector& regionStart,
                                               const Rcpp::IntegerVector& regionEnd,
                                               const Rcpp::CharacterVector& seqstring,
                                               const Rcpp::IntegerVector& seqStart,
                                               const Rcpp::IntegerVector& seqEnd,
                                               const Rcpp::LogicalVector& remove_nonunique,
                                               const Rcpp::IntegerVector& clip_until_nbg,
                                               const Rcpp::NumericVector& max_protect_frac,
                                               const Rcpp::NumericVector& max_bisC_meth,
                                               const Rcpp::IntegerVector& min_bisC_size,
                                               const Rcpp::IntegerVector& mapqMin,
                                               const Rcpp::IntegerVector& mapqMax,
                                               const Rcpp::CharacterVector& alignerUsed,
                                               const Rcpp::CharacterVector& SMFenzymeUsed,
                                               const Rcpp::IntegerVector& min_frag_data_len,
                                               const Rcpp::NumericVector& min_frag_data_dens,
                                               const Rcpp::LogicalVector& data_as_rle);

// [[Rcpp::export]]
Rcpp::List fetch_cooc_ctable_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::IntegerVector& max_spac,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax,
                                           const Rcpp::CharacterVector& alignerUsed,
                                           const Rcpp::CharacterVector& SMFenzymeUsed,
                                           const Rcpp::IntegerVector& min_frag_data_len,
                                           const Rcpp::NumericVector& min_frag_data_dens);


// [[Rcpp::export]]
Rcpp::List fetch_data_matrix_from_bams_cpp(const Rcpp::CharacterVector& whichContext,
                                           const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax,
                                           const Rcpp::CharacterVector& alignerUsed,
                                           const Rcpp::CharacterVector& SMFenzymeUsed,
                                           const Rcpp::IntegerVector& min_frag_data_len,
                                           const Rcpp::NumericVector& min_frag_data_dens);

// [[Rcpp::export]]
Rcpp::List fetch_protect_stats_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                             const Rcpp::CharacterVector& regionChr,
                                             const Rcpp::IntegerVector& regionStart,
                                             const Rcpp::IntegerVector& regionEnd,
                                             const Rcpp::CharacterVector& seqstring,
                                             const Rcpp::IntegerVector& seqStart,
                                             const Rcpp::IntegerVector& seqEnd,
                                             const Rcpp::LogicalVector& remove_nonunique,
                                             const Rcpp::IntegerVector& clip_until_nbg,
                                             const Rcpp::NumericVector& max_protect_frac,
                                             const Rcpp::NumericVector& max_bisC_meth,
                                             const Rcpp::IntegerVector& min_bisC_size,
                                             const Rcpp::IntegerVector& mapqMin,
                                             const Rcpp::IntegerVector& mapqMax,
                                             const Rcpp::CharacterVector& alignerUsed,
                                             const Rcpp::CharacterVector& SMFenzymeUsed,
                                             const Rcpp::IntegerVector& min_frag_data_len,
                                             const Rcpp::NumericVector& min_frag_data_dens);

#endif