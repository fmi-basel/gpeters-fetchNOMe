#ifndef _region_data_
#define _region_data_

#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

#include "fragData-class.h"
#include "utils_globvars.h"
#include "protectStats-class.h"

using namespace std;


// container for storing all fragments overlapping a region of interest
class regionData{
	string regChr;
	int regStart;
	int regEnd;
	
	unordered_map<string, fragData > fDataMap; // map: qname (string) -> fragment (class fragData)
	
public:
  
  /////  constructors/destructors ///// 
	regionData();
	regionData(const string& regchr,
            const int& regstart,
            const int& regend);
	~regionData();
	
	/////  functions for adding new data ///// 
	bool addFragContextData(const string& bampref,
                         const string& qname,
                         const int& fragstart,
                         const int& fragend,
                         const string& fragconfig,
                         const vector<int >& rposvec,
                         const vector<bool >& protectvec,
                         const vector<bool >& isonplusstrandvec,
                         const bool& isSecond,// data for same positions on the second read are ignored
                         fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
	);
	
	void print();
	
	/////  functions for collecting statistics ///// 
	void addCountsToCoocTable(coocCntTable& cTable);
	
	
	void addCountsToProtectStatsTable(protectStatsTable& pStatsTable);
	
	
	///// functions for post-processing data /////
	// clip fragments to get rid of partial footprints
	int clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                     const double& max_protect_frac); // keeps fragment if percentage of 1s in a fragment is above this number)
	
	
	
	// get set of qnames whose meth profile is not unique for removal
	set<string > getNonUniqueQnames();
	
	// get set of qnames with failed bisulfite conversion
	set<string > getFailedBisConvQnames(const double& max_bisC_meth,
                                     const int& min_bisC_size);
	
	// remove fragments from the map. returns new size of the map
	int eraseQnames(const set<string >& qnames);
	
	// remove fragments which contain NO data. used after clipping
	int eraseEmptyFrags(fragMapType whichMap);
	
	// get set of qnames whose data length or data density is too short
	set<string > getInsuffDataQnames(const int& min_frag_data_len,
                                  const double& min_frag_data_dens,
                                  const fragMapType& whichMap);
	
	
	vector<int > postProcessFragData(const bool& remove_nonunique,
                                  const int& clip_until_nbg,
                                  const double& max_protect_frac,
                                  const double& max_bisC_meth,
                                  const int& min_bisC_size,
                                  const int& min_frag_data_len,
                                  const double& min_frag_data_dens,
                                  const fragMapType& whichMap);
	
	
	
	int size();
	
	
	// convert data in region into a Matrix and return
	Rcpp::IntegerMatrix getDataMatrix(fragMapType whichMap);
	
	// get List with fetched reads/fragments for R 
	Rcpp::List getFragDataList(fragMapType whichMap);
	
	// get List with fetched reads/fragments for R with Rle compressed data vectors
	Rcpp::List getFragDataRle(fragMapType whichMap);
	
};



#endif