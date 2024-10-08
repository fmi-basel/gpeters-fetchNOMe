#ifndef _fetch_regdata_from_bam_hpp_
#define _fetch_regdata_from_bam_hpp_

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

#ifdef __cplusplus
extern "C" {
#endif
#include "utilities.h"
#ifdef __cplusplus
}
#endif

using namespace std;


typedef struct { // for use in callback function. contains pointers to required container objects
	refSeqInfo *refseq_info; // pointer to an object with information about reference sequence
	regionData *reg_data; // point to an object which stores data in region
	int mapqMin; // minimum mapping quality score
	int mapqMax; // maximum mapping quality score
	int absIsizeMin; // minimum absolute tlen for paired reads
	int absIsizeMax; // maximum absolute tlen for paired reads
	int min_read_size; // minimum read length
	int max_read_size; // maximum read length
	string prefix; // prefix for qnames to avoid problems with identical qnames in different bams
} obj_pnts;


// fetch data from bam file
bool fetch_data_from_bam(const string& bamfile,
                         const string& regionChr,
                         const int& regionStart,
                         const int& regionEnd,
                         const alignerType& aligner,
                         obj_pnts* pnts);

// bam_fetch callback function for bam_fetch for data retrieval from BAM files produced by QuasR and Bismark
static int collectRegionData_QsBism(const bam1_t *hit, void *data);

// bam_fetch callback function for bam_fetch for data retrieval from BAM files produced by BISCUIT
static int collectRegionData_Bisq(const bam1_t *hit, void *data);



// checks validity of an alignment
bool isValidAln(const bam1_t *hit,const obj_pnts* pdata);

// print alignment. for debugging
void print_bs_alignment(const bam1_t *hit,const obj_pnts* pdata);

#endif