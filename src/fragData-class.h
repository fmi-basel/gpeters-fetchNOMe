#ifndef _frag_data_
#define _frag_data_


#include <Rcpp.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cstdint>
#include <stdbool.h>
#include <stdlib.h>
#include "refSeqInfo-class.h"
#include "utils_globvars.h"
#include "protectStats-class.h"
#include "coocCntTable-class.h""


using namespace std;

// container class for protection data based on GCH methylation, endogenous methylation at WCG
// methylation for estimating bisulfite conversion efficiency and all other C methylation

class fragData{
	// keys in maps are 1-based positions on reference
	
	// fragment 1-based coordinates on reference
	int fragStart;
	int fragEnd;
	
	// configuration of a fragment, i.e. R1+R2- etc.
	string fragConfig;
	// file prefix where the fragment came from
	string prefix;
	
	// protection at GCH sites
	map<int, bool > gchProtect; // map: 1-based position -> protection {0 - unprotected, 1 - protected}
	map<int, bool > gchIsOnPlusStrand; // 1 - "+" strand, 0 - "-" strand.
	// The reason that the strand is encoded as an vector is that fragData is a container for a fragment consisting of both reads R1 and R2
	
	// methylation of endogenous CpG at WCG sites
	map<int , bool > wcgMeth;
	map<int , bool > wcgIsOnPlusStrand;
	
	// methylation at sites for estimating bisulfite conversion efficiency
	map<int , bool > bisCMeth;
	map<int , bool > bisCIsOnPlusStrand;
	
	// methylation of other C's
	map<int , bool > otherCMeth;
	map<int , bool > otherCIsOnPlusStrand;
public:
	
	
	fragData();
	
	// constructors/destructors
	fragData(const string& bampref,
          const int& fragstart,
          const int& fragend,
          const string& fragconfig,
          const vector<int >& gch_pos,
          const vector<bool >& gch_protect,
          const vector<bool >& gch_onplusstrand,
          const vector<int >& wcg_pos,
          const vector<bool >& wcg_protect,
          const vector<bool >& wcg_onplusstrand,
          const vector<int >& bisc_pos,
          const vector<bool >& bisc_protect,
          const vector<bool >& bisc_onplusstrand,
          const vector<int >& otherc_pos,
          const vector<bool >& otherc_protect,
          const vector<bool >& otherc_onplusstrand
	);
	
	~fragData();
	
	// get vector of positions in a map
	vector<int > getPositions(fragMapType whichMap);
	
	
	// get data at position
	bool getDataAt(int rpos,
                fragMapType whichMap);
	
	// get isOnPlus
	bool isDataOnPlusStrand(int rpos,
                         fragMapType whichMap);
	
	
	
	// add data into the map 
	bool addData(const string& bampref,
              const int& fragstart,
              const int& fragend,
              const string& fragconfig,
              const vector<int >& rposvec,
              const vector<bool >& protectvec,
              const vector<bool >& isonplusstrandvec,
              const bool& isSecond,// data for same positions on the second read are ignored
              fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
	);
	
	// convert data to strings
	vector<string > getStringData();
	
	string getFragConfig();
	
	// // get most left position of the first point with data
	// int getFragDataStart(fragMapType whichMap);
	// 
	// // get most right position of the last point with data
	// int getFragDataEnd(fragMapType whichMap);
	
	// get positions where the data start, end and number of informative positions
	std::vector<int > getFragDataStartEndNpoints(fragMapType whichMap);
	
	
	
	
	// get the data, npoints and left and right datapoints;
	Rcpp::List getFragDataVector(fragMapType whichMap);
	
	
	// get the Rle compressed data, npoints and left and right datapoints;
	Rcpp::List getFragDataRle(fragMapType whichMap);
	
	
	// print
	void print();
	
	// 3. function to add cooccurrence counts into provided table
	void addCountsToCoocTable(coocCntTable& cTable);
	
	
	// function to add counts to protection statistics table
	void addCountsToProtectStatsTable(protectStatsTable& pStatsTable);
	
	
	
	// 2. function to clip data. returns 1 if fragment was clipped or 0 otherwise
	int clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                     const double& max_protect_frac // keeps fragment if percentage of 1s in a fragment is above this number
	);
	
	
	double get_bisC_meth(const int& min_size);
	
	bool isEmpty(fragMapType whichMap);
	
	// get string representation of fragData
	string getStringEncodingForUniqueness();
	
};




#endif