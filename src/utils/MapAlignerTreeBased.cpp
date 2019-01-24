#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SpanningGraph.h>
#include <boost/regex.hpp>
#include <algorithm>
#include <numeric>
#include <math.h>  

using namespace OpenMS;
using namespace std;


class MapAlignerTreeBased:
  public TOPPMapAlignerBase
{
public:
  MapAlignerTreeBased():
	TOPPMapAlignerBase("MapAlignerTreeBased","Tree-based alignment with identification algorithm.")
	{
	}

private:
	void registerOptionsAndFlags_()
	{
	  String file_formats = "featureXML,consensusXML";
	  registerInputFileList_("in", "<files>", StringList(), "Input files to align (all must have the same file type)", true);
         setValidFormats_("in", ListUtils::create<String>(file_formats));
         registerOutputFile_("out", "<file>", "", "Output file.", false);
         setValidFormats_("out", ListUtils::create<String>(file_formats));
         registerSubsection_("algorithm", "Algorithm parameters section");
         registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
         registerFlag_("keep_subelements", "For consensusXML input only: If set, the sub-features of the inputs are transferred to the output.");
      } 
  
  ExitCodes main_(int, const char**)
  {
    StringList input_files = getStringList_("in");
    String output_file = getStringOption_("out");
    vector<ConsensusMap> maps;
    vector<FeatureMap> fmaps;
    FileTypes::Type in_type = FileHandler::getType(input_files[0]); 
    
    if (in_type == FileTypes::FEATUREXML)
    {
	  for(unsigned int  i = 0; i < input_files.size(); i++)
	  {
	    FeatureMap f;
	    ConsensusMap c;
           FeatureXMLFile().load(input_files[i],f);
           //fmaps.push_back(f);
	    MapConversion::convert(0, f, c, -1);
	    c.getColumnHeaders()[0].filename = input_files[i]; //get filenames
           c.applyMemberFunction(&UniqueIdInterface::setUniqueId);
           c.getColumnHeaders()[0].size = c.size();
	    maps.push_back(c);   
	  }
    }
    
    else if (in_type == FileTypes::CONSENSUSXML)
    {
      int index = 0;
      for(StringList::const_iterator it = input_files.begin(); it!=input_files.end(); ++it)
      {
        ConsensusMap c;
        ConsensusXMLFile().load(*it,c);
        c.getColumnHeaders()[index].filename = input_files[index]; //get filenames
        c.getColumnHeaders()[index].size = maps[index].size();
        c.applyMemberFunction(&UniqueIdInterface::setUniqueId);
        maps.push_back(c);
        index++;
      }	
    }
        
    vector<unsigned int> ori_index;
    vector<pair<int,String>> files;
        
    vector<vector<double>> M;
    ConsensusMap out;
    
    for (unsigned int i = 0; i < maps.size(); i++)
    {
      vector<double> row;
      for (unsigned int j = 0; j < maps.size(); j++)
      {
        if (i==j) {row.push_back(DBL_MAX);}
        else {row.push_back(0);}
      }
      M.push_back(row);
    }
    
    computeMetric(M,maps);   
    vector<pair<pair<int,int>,float>> queue;
    computeSpanningTree(M,queue); 
    alignSpanningTree(queue,maps,input_files,out); 
    
    ConsensusXMLFile().store(output_file, out); 

    return EXECUTION_OK;
  
  }
  
  
  
  Param getSubsectionDefaults_(const String & section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return TOPPMapAlignerBase::getModelDefaults("b_spline");
    }

    return Param(); 
  }
  
  

  //Fill in matrix
  void computeMetric(vector<vector<double>>& matrix, vector<ConsensusMap>& maps) const
  { 
  
  vector<map<pair<AASequence,int>,DoubleList>> all_seq;
  
  for (Size m = 0; m < maps.size(); m++) 
  {
  
    map<pair<AASequence,int>,DoubleList> seq_rts;
    vector<PeptideIdentification> un_pep = maps[m].getUnassignedPeptideIdentifications();
    for (vector<PeptideIdentification>::iterator pep_it = un_pep.begin(); pep_it!=un_pep.end();pep_it++)
    {
    
      pair<AASequence,int> seq_ch = make_pair(pep_it->getHits()[0].getSequence(),pep_it->getHits()[0].getCharge());
      seq_rts[seq_ch].push_back(pep_it->getRT());
    }
    
    for (vector<ConsensusFeature>::iterator c_it = maps[m].begin(); c_it!=maps[m].end();c_it++) 
    {
      
      if (!c_it->getPeptideIdentifications().empty())
      { 
        
        for (vector<PeptideIdentification>::iterator p_it = 
             c_it->getPeptideIdentifications().begin(); 
             p_it!=c_it->getPeptideIdentifications().end();++p_it)
          {
            if (!p_it->getHits().empty())
            {
              //Writes peptide hit with the highest score
              p_it->sort();
              pair<AASequence,int> seq_ch = make_pair(p_it->getHits()[0].getSequence(),p_it->getHits()[0].getCharge());
              seq_rts[seq_ch].push_back(p_it->getRT());
              
            }
          }
       }
    }
    all_seq.push_back(seq_rts);
  
  }
  
  for(unsigned int i = 0; i < all_seq.size()-1; i++)
  {
    for(unsigned int j = i; j < all_seq.size(); j++) 
    {
      vector<float> map1;
      vector<float> map2;
      
      for (map<pair<AASequence,int>,DoubleList>::iterator a_it = all_seq[i].begin(); a_it != all_seq[i].end(); ++a_it)
      {
        for (map<pair<AASequence,int>,DoubleList>::iterator b_it = all_seq[j].begin(); b_it != all_seq[j].end(); ++b_it)
        {
          if(a_it->first==b_it->first)
          {
            pair<double,double> min = closestMatch(a_it->second,b_it->second);
            map1.push_back(min.first);
            map2.push_back(min.second);
          }
        }
      }
      if (map1.size()>2)
      {
        
        double pearson = Math::pearsonCorrelationCoefficient(map1.begin(), map1.end(), map2.begin(), map2.end());
        LOG_INFO << "Found " << map1.size() << " matching peptides for " << i << " and " << j << endl;
        //Math::computeRank(map1);
        //Math::computeRank(map2);        
        //double rpearson = Math::rankCorrelationCoefficient(map1.begin(), map1.end(), map2.begin(), map2.end());
        if (!isnan(pearson))
        {
          matrix[i][j] = 1-abs(pearson); 
          matrix[j][i] = 1-abs(pearson);
          LOG_INFO << matrix[i][j] << endl;
        }
        else
        {
          matrix[i][j] = 1;
          matrix[j][i] = 1;
          LOG_INFO << matrix[i][j] << endl;
        }
      }
      else 
      {
        matrix[i][j] = 2; //no correlation
        matrix[j][i] = 2; 
        LOG_INFO << matrix[i][j] << endl;
      } 
    }
  }
  
  
}


//find closest match of same peptides in maps
pair<double,double> closestMatch(DoubleList a, DoubleList b) const
{ 
  double RT1 = DBL_MAX;
  double RT2 = DBL_MAX;
  double min_dist = DBL_MAX;
  for (DoubleList::iterator a_it = a.begin(); a_it != a.end(); ++a_it )
  {
    for(DoubleList::iterator b_it = b.begin(); b_it != b.end(); ++b_it)
    {
      double dist = abs(*a_it - *b_it);
      if(dist<min_dist) {RT1 = *a_it; RT2 = *b_it;}
    }
  
  }
  return make_pair(RT1,RT2);

}


//compute MST for the tree-based alignment
void computeSpanningTree(vector<vector<double>> matrix, 
                         vector<pair<pair<int,int>,float>>& queue)
                         
{ 
  vector<pair<pair<int,int>,float>> sortedEdges;
  for (unsigned int i = 0; i < matrix.size()-1; i++)
  {
    for (unsigned int j = i+1; j < matrix.size(); j++)
    {
      sortedEdges.push_back(make_pair(make_pair(i,j),matrix[i][j]));
    }
  }
  
  std::sort(sortedEdges.begin(), sortedEdges.end(), sortByScore);
 
  int size = matrix.size();
  queue.push_back(sortedEdges[0]);
  sortedEdges.erase(sortedEdges.begin());
  
  while(!sortedEdges.empty())
  {
    SpanningGraph tree(size);
    for (unsigned int i = 0; i < queue.size(); i++)
    {
      tree.addEdge(queue[i].first.first,queue[i].first.second);
    }
    tree.addEdge(sortedEdges[0].first.first,sortedEdges[0].first.second);
    if (tree.containsCycle()) {}
    else  {queue.push_back(sortedEdges[0]);}
    sortedEdges.erase(sortedEdges.begin());
  }
  
  SpanningGraph aligner(size);
  LOG_INFO << "Minimum spanning graph: " << endl;
  for (unsigned int i = 0; i < queue.size(); i++)
  {
    LOG_INFO << queue[i].first.first << " <-> " << queue[i].first.second << ", weight: " << queue[i].second << endl;
    aligner.addEdge(queue[i].first.first,queue[i].first.second);
  }

 
}

//sorting function for queue
static bool sortByScore(const pair<pair<int,int>,float> &lhs, const pair<pair<int,int>,float> &rhs) { return lhs.second < rhs.second; }

//alignment util
void align(vector<ConsensusMap>& to_align, vector<TransformationDescription>& transformations)
{
  
  MapAlignmentAlgorithmIdentification algorithm;
  Param algo_params = getParam_().copy("algorithm:", true);
  algorithm.setParameters(algo_params);
  algorithm.setLogType(log_type_);
  
  algorithm.align(to_align,transformations);
  Param model_params = getParam_().copy("model:", true);
  String model_type = model_params.getValue("type");
  if (model_type != "none")
  {
      model_params = model_params.copy(model_type + ":", true);
      for (vector<TransformationDescription>::iterator it = 
             transformations.begin(); it != transformations.end(); ++it)
      {	
        it->fitModel(model_type, model_params);
      }
  }
  
  for (unsigned int i = 0; i < transformations.size(); i++)
  {
      MapAlignmentTransformer::transformRetentionTimes(to_align[i], transformations[i]);
      addDataProcessing_(to_align[i], getProcessingInfo_(DataProcessing::ALIGNMENT));
  }
 
}

//Main alignment function
void alignSpanningTree(vector<pair<pair<int,int>,float>>& queue, vector<ConsensusMap>& maps,
                       StringList input_files, ConsensusMap& out_map)
{  
 
  for (unsigned int i = 0; i < queue.size(); i++)
  {
    vector<ConsensusMap> to_align;
    int A = queue[i].first.first;
    int B = queue[i].first.second;
    to_align.push_back(maps[A]);
    to_align.push_back(maps[B]);
    
    vector<TransformationDescription> transformations;
    for (unsigned int j = 0; j < to_align.size(); j++)
    {
      TransformationDescription t;
      transformations.push_back(t); 
    }
    align(to_align,transformations);
  
    //Grouping step
  
    ConsensusMap out;
     
    StringList ms_run_locations;
    for (Size j = 0; j < to_align.size(); ++j)
    {
        to_align[j].updateRanges();
        StringList ms_runs;
        to_align[j].getPrimaryMSRunPath(ms_runs);
        ms_run_locations.insert(ms_run_locations.end(), ms_runs.begin(), ms_runs.end());
    }
 
    FeatureGroupingAlgorithmQT grouping;
    
    out.getColumnHeaders()[0].filename = input_files[A];
    out.getColumnHeaders()[0].size = maps[A].size();
    out.getColumnHeaders()[0].unique_id = maps[A].getColumnHeaders()[0].unique_id;
    out.getColumnHeaders()[1].filename = input_files[B];
    out.getColumnHeaders()[1].size = maps[B].size();
    out.getColumnHeaders()[1].unique_id = maps[B].getColumnHeaders()[1].unique_id;
    
    grouping.group(to_align,out);
    
    grouping.transferSubelements(to_align,out);
  
    out.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    addDataProcessing_(out,getProcessingInfo_(DataProcessing::FEATURE_GROUPING));

    out.sortPeptideIdentificationsByMapIndex();
  
    
    if (i < queue.size()-1)
    {
      for (unsigned int q = i+1; q < queue.size(); q++)
      {
        if (queue[q].first.first==B) {queue[q].first.first=A;}
        else if (queue[q].first.second==B) {queue[q].first.second=A;}
      }
    }
     
     
    maps[A] = out;
    maps[B] = out;
        
    map<Size, UInt> num_consfeat_of_size;
    for (ConsensusMap::const_iterator cmit = out.begin();
         cmit !=out.end(); ++cmit)
    {
      ++num_consfeat_of_size[cmit->size()];
    }

    LOG_INFO << "Number of consensus features:" << endl;
    for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin();
         i != num_consfeat_of_size.rend(); ++i)
    {
      LOG_INFO << "  of size " << setw(2) << i->first << ": " << setw(6) 
               << i->second << endl;
    }
    LOG_INFO << "  total:      " << setw(6) << out.size() << endl;
    
    to_align.clear();
    out_map = out;
    out.clear();
  }

} 

};


int main(int argc, const char** argv)
{
  MapAlignerTreeBased tool;
  return tool.main(argc,argv);
}
