#include "regionData-class.h"

regionData::regionData(){};
regionData::regionData(const string& regchr,
           const int& regstart,
           const int& regend){
  regChr = regchr;
  regStart = regstart;
  regEnd = regend;
};
regionData::~regionData(){};

bool regionData::addFragContextData(const string& bampref,
                        const string& qname,
                        const int& fragstart,
                        const int& fragend,
                        const string& fragconfig,
                        const vector<int >& rposvec,
                        const vector<bool >& protectvec,
                        const vector<bool >& isonplusstrandvec,
                        const bool& isSecond,// data for same positions on the second read are ignored
                        fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
){
  // check if qname exists.
  // no - create new fragment, add data and insert into fDataMap
  // yes - add additional data into corresponding fragment
  
  if(fDataMap.find(qname) == fDataMap.end()){ // if not found a query name then insert a new
    fragData frg;
    frg.addData(bampref,
                fragstart,
                fragend,
                fragconfig,
                rposvec,
                protectvec,
                isonplusstrandvec,
                isSecond,// data for same positions on the second read are ignored
                whichMap);
    fDataMap.insert(make_pair(qname, frg));
    
  }
  else { // if found then add to existing fragment.
    fDataMap[qname].addData(bampref,
                            fragstart,
                            fragend,
                            fragconfig,
                            rposvec,
                            protectvec,
                            isonplusstrandvec,
                            isSecond,// data for same positions on the second read are ignored
                            whichMap);
  }
  return(1);
};

void regionData::print(){
  Rcpp::Rcout<<"=========== Data in "<<regChr<<":"<<regStart<<"-"<<regEnd<<" ==========="<<endl;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    Rcpp::Rcout<<"==== fragment name: "<<it->first<<" ===="<<endl;
    (it->second).print();
  }
  
  Rcpp::Rcout<<"========================================================================"<<endl;
};


void regionData::addCountsToCoocTable(coocCntTable& cTable){
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    (it->second).addCountsToCoocTable(cTable);
  }
};


void regionData::addCountsToProtectStatsTable(protectStatsTable& pStatsTable){
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    (it->second).addCountsToProtectStatsTable(pStatsTable);
  }
  
}


// clip fragments to get rid of partial footprints
int regionData::clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                    const double& max_protect_frac) // keeps fragment if percentage of 1s in a fragment is above this number)
{
  int nFragsClipped = 0;
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    nFragsClipped += (it->second).clipProtectData(clip_until_nbg,max_protect_frac);
  }
  
  return(nFragsClipped);
};


// get set of qnames whose meth profile is not unique for removal
set<string > regionData::getNonUniqueQnames(){
  set<string > nUqnames;
  unordered_set<string > fragdatenc;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    
    string enc = (it->second).getStringEncodingForUniqueness();
    // for testing
    // (it->second).print();
    // Rcpp::Rcout<<enc<<endl;
    
    if(fragdatenc.find(enc) != fragdatenc.end()){ // encoding exists
      nUqnames.insert(it->first); // add qname to set of non-unique fragments
    } else{ // new encoding
      fragdatenc.insert(enc); 
    }
  }
  return(nUqnames);
  
};

// get set of qnames with failed bisulfite conversion
set<string > regionData::getFailedBisConvQnames(const double& max_bisC_meth,const int& min_bisC_size){
  set<string > failed_qnames;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    double biscmeth = (it->second).get_bisC_meth(min_bisC_size);
    if(biscmeth > max_bisC_meth)
      failed_qnames.insert(it->first);
  }
  return(failed_qnames);
  
};

// remove fragments from the map. returns new size of the map
int regionData::eraseQnames(const set<string >& qnames){
  for(set<string >::iterator it=qnames.begin(); it != qnames.end(); ++it)
    fDataMap.erase(*it);
  return(fDataMap.size());
};


int regionData::eraseEmptyFrags(fragMapType whichMap){
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();){
    if((it->second).isEmpty(whichMap)){
      it = fDataMap.erase(it);
    } else{
      ++it;
    }
  }
  return(fDataMap.size());
};

set<string > regionData::getInsuffDataQnames(const int& min_frag_data_len,
                                             const double& min_frag_data_dens,
                                             const fragMapType& whichMap){
  set<string > failed_qnames;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    std::vector<int > dataPosNpoints = (it->second).getFragDataStartEndNpoints(whichMap);
    
    int fragDataStart = dataPosNpoints[0];
    int fragDataEnd = dataPosNpoints[1];
    int npoints = dataPosNpoints[2];
    int fragDataLen = fragDataEnd - fragDataStart + 1;
    
    if(fragDataLen == 0 || npoints == 0){
      failed_qnames.insert(it->first);
      continue;
    }
    
    double fragDataDens = static_cast<double>(npoints) / fragDataLen;
    
    if(fragDataLen < min_frag_data_len || fragDataDens < min_frag_data_dens)
      failed_qnames.insert(it->first);
  }
  return(failed_qnames);
};


vector<int > regionData::postProcessFragData(const bool& remove_nonunique, // remove non-unique?
                                             const int& clip_until_nbg, // clipping partial footprints at the edges
                                             const double& max_protect_frac, //
                                             const double& max_bisC_meth,
                                             const int& min_bisC_size,
                                             const int& min_frag_data_len,
                                             const double& min_frag_data_dens,
                                             const fragMapType& whichMap){
  ///// 1. filtering non-unique fragments /////
  // get qnames on non-unique fragments
  set<string > nonUnfrags;
  if(remove_nonunique)
    nonUnfrags = getNonUniqueQnames();
  int n_nonUnique = nonUnfrags.size();
  
  
  ///// 2. filtering fragments with failed bisulfite conversion, i.e. too high bisC methylation /////
  // get qnames with failed bisulfite conversion
  set<string > failedBisConvfrags;
  if(max_bisC_meth < 1.0)
    failedBisConvfrags = getFailedBisConvQnames(max_bisC_meth,min_bisC_size);
  int n_bisfailed = failedBisConvfrags.size();
  
  // erase fragments which are nonunique or failed bisConv
  set<string > remove_qnames(nonUnfrags);
  remove_qnames.insert(failedBisConvfrags.begin(), failedBisConvfrags.end());
  eraseQnames(remove_qnames);
  // empty the set with qnames
  remove_qnames.erase(remove_qnames.begin(),remove_qnames.end());
  
  ///// 3. clip fragments /////
  int n_clipped = 0;
  if(clip_until_nbg > 0){
    n_clipped = clipProtectData(clip_until_nbg,max_protect_frac);
    // // remove fragments which were erased by clipping
    // eraseEmptyFrags(fragMapType whichMap);
  }
  
  //// 4. filtering fragments by data length, i.e. genomic distance from most left and most right informative (containing data) positions
  set<string > insuffDatafrags = getInsuffDataQnames(min_frag_data_len,
                                                     min_frag_data_dens,
                                                     whichMap);
  int n_insuffFailed = insuffDatafrags.size();
  
  // erase fragments which contain insufficient data
  int n_analysed = eraseQnames(insuffDatafrags);
  
  std::vector<int > vec_out={n_nonUnique,n_bisfailed,n_clipped,n_insuffFailed,n_analysed};
  return(vec_out);
};




int regionData::size(){
  return(fDataMap.size());
};


// convert data in region into a Matrix and return
Rcpp::IntegerMatrix regionData::getDataMatrix(fragMapType whichMap){
  
  // create Matrix filled with NA_Integer
  Rcpp::IntegerMatrix matr(fDataMap.size(),regEnd - regStart);
  fill(matr.begin(),matr.end(),NA_INTEGER);
  // set column names
  Rcpp::CharacterVector colnm(regEnd - regStart);
  for(int pos=regStart + 1; pos <= regEnd; ++pos)
    colnm[pos - regStart - 1] = to_string(pos);
  Rcpp::colnames(matr) = colnm;
  // set rownames
  vector<string > rownm;
  for(unordered_map<string, fragData >::iterator it=fDataMap.begin(); it != fDataMap.end(); ++it)
    rownm.push_back(it->first);
  Rcpp::rownames(matr) = Rcpp::wrap(rownm);
  
  if(whichMap != allC){
    for(int rowi=0; rowi < rownm.size();++rowi){
      // get positions in a fragment
      vector<int > pos_vec = fDataMap[rownm[rowi]].getPositions(whichMap);
      for(int posi=0; posi < pos_vec.size(); ++posi){
        if(pos_vec[posi] >= regStart && pos_vec[posi] < regEnd)
          matr(rowi,pos_vec[posi] - regStart) = fDataMap[rownm[rowi]].getDataAt(pos_vec[posi],whichMap);
      }
    }
  } else{
    // if data for all C's is requested then loop across all contexts and fill the matrix
    vector<fragMapType > frag_map_vec = {GCH, WCG, bisC, otherC};
    for(int rowi=0; rowi < rownm.size();++rowi){
      // loop across context map type
      for(int mtpi = 0; mtpi < frag_map_vec.size(); ++mtpi){
        // get positions in a fragment
        vector<int > pos_vec = fDataMap[rownm[rowi]].getPositions(frag_map_vec[mtpi]);
        for(int posi=0; posi < pos_vec.size(); ++posi){
          if(pos_vec[posi] >= regStart && pos_vec[posi] < regEnd){
            // get data at position
            int posdat = fDataMap[rownm[rowi]].getDataAt(pos_vec[posi],frag_map_vec[mtpi]);
            // assign to matrix element, HOWEVER invert for context GCH as we store protection values which are 0-unprotected->methylated and 1-protected->unmethylated
            matr(rowi,pos_vec[posi] - regStart) = frag_map_vec[mtpi] == GCH ? 1-posdat : posdat;
          }
          
        }
      }
    }
  }
  
  return(matr);
};


Rcpp::List regionData::getFragDataList(fragMapType whichMap){
  // start (genomic position of the first data point. note, not the alignment start!)
  // end last position of the data points
  // fragConfigs configurations of fragments
  // qname - name of the fragment
  // fragData_list list which contains vectors with 0 (accessible), 1(protected) and NA, non-informative
  
  
  Rcpp::IntegerVector fragStarts; // start positions are 1-based on reference
  Rcpp::IntegerVector fragEnds; // start positions are 1-based on reference
  
  Rcpp::IntegerVector nDataPoints; // number of informative points
  Rcpp::IntegerVector startDataPoints; // data points at start positions
  Rcpp::IntegerVector endDataPoints; // data points at end positions
  
  Rcpp::CharacterVector fragConfigs; // configurations of configs
  
  Rcpp::CharacterVector fragQnames; // fragment names
  
  Rcpp::List fragData_list;
  
  
  for(unordered_map<string, fragData >::iterator it=fDataMap.begin(); it != fDataMap.end(); ++it){
    // get number of informative points in a fragment
    Rcpp::List fragdat = (it->second).getFragDataVector(whichMap);
    int npoints = Rcpp::as<int >(fragdat["nDataPoints"]);
    
    int fragDataStart = Rcpp::as<int >(fragdat["fragDataStart"]);
    int fragDataEnd = Rcpp::as<int >(fragdat["fragDataEnd"]);
    
    int fragDataLen = fragDataEnd - fragDataStart + 1;
    
    fragQnames.push_back(it->first);
    fragStarts.push_back(fragDataStart + 1); // +1 is to convert to 1-based positions
    fragEnds.push_back(fragDataEnd + 1);
    fragConfigs.push_back((it->second).getFragConfig());
    nDataPoints.push_back(fragdat["nDataPoints"]);
    startDataPoints.push_back(fragdat["startData"]);
    endDataPoints.push_back(fragdat["endData"]);
    fragData_list.push_back(fragdat["fragDataVec"]);
  }
  
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("start") = fragStarts,
                                           Rcpp::Named("end") = fragEnds,
                                           Rcpp::Named("fragConfigs") = fragConfigs,
                                           Rcpp::Named("qnames") = fragQnames,
                                           Rcpp::Named("nDataPoints") = nDataPoints,
                                           Rcpp::Named("startDataPoint") = startDataPoints,
                                           Rcpp::Named("endDataPoint") = endDataPoints,
                                           Rcpp::Named("fragData_list") = fragData_list);
  return(data_out);
};


Rcpp::List regionData::getFragDataRle(fragMapType whichMap){
  // start (genomic position of the first data point. note, not the alignment start!)
  // end last position of the data points
  // fragConfigs configurations of fragments
  // qname - name of the fragment
  // fragData_list list which contains Rle compressed vectors with 0 (accessible), 1(protected) and NA, non-informative
  
  
  Rcpp::IntegerVector fragStarts; // start positions are 1-based on reference
  Rcpp::IntegerVector fragEnds; // start positions are 1-based on reference
  
  Rcpp::IntegerVector nDataPoints; // number of informative points
  Rcpp::IntegerVector startDataPoints; // data points at start positions
  Rcpp::IntegerVector endDataPoints; // data points at end positions
  
  Rcpp::CharacterVector fragConfigs; // configurations of configs
  
  Rcpp::CharacterVector fragQnames; // fragment names
  
  Rcpp::List fragData_list;
  
  
  
  for(unordered_map<string, fragData >::iterator it=fDataMap.begin(); it != fDataMap.end(); ++it){
    
    Rcpp::List fragdat = (it->second).getFragDataRle(whichMap);
    int npoints = Rcpp::as<int >(fragdat["nDataPoints"]);
    
    int fragDataStart = Rcpp::as<int >(fragdat["fragDataStart"]);
    int fragDataEnd = Rcpp::as<int >(fragdat["fragDataEnd"]);
    
    int fragDataLen = fragDataEnd - fragDataStart + 1;
    
    fragQnames.push_back(it->first);
    fragStarts.push_back(fragDataStart + 1); // +1 is to convert to 1-based positions
    fragEnds.push_back(fragDataEnd + 1);
    fragConfigs.push_back((it->second).getFragConfig());
    nDataPoints.push_back(fragdat["nDataPoints"]);
    startDataPoints.push_back(fragdat["startData"]);
    endDataPoints.push_back(fragdat["endData"]);
    
    fragData_list.push_back(fragdat["rle_data"]);
  }
  
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("start") = fragStarts,
                                           Rcpp::Named("end") = fragEnds,
                                           Rcpp::Named("fragConfigs") = fragConfigs,
                                           Rcpp::Named("qnames") = fragQnames,
                                           Rcpp::Named("nDataPoints") = nDataPoints,
                                           Rcpp::Named("startDataPoint") = startDataPoints,
                                           Rcpp::Named("endDataPoint") = endDataPoints,
                                           Rcpp::Named("fragData_list") = fragData_list);
  return(data_out);
};


