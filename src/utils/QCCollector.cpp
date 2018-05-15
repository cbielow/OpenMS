// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Anton Haberland, Leo Wurth
// $Authors: Anton Haberland, Leo Wurth
// --------------------------------------------------------------------------
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetrics.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;


class OPENMS_DLLAPI QCCollector:
  public TOPPBase
{
public:
  QCCollector():
	TOPPBase("QCCollector","Will collect several mzTabs from several utils.",false)
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
	  registerInputFileList_("in_protein_quantifier_peptide","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_protein_quantifier_protein","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_IDMapper","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_map_RT_transformer","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_internal_calibration","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_rawfiles_false_discovery_rate","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_post_false_discovery_rate","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_feature_linker_unlabeledQT","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_contaminant_database","<files>",StringList(), "Input files",false,false);
    setValidFormats_("in_protein_quantifier_peptide", ListUtils::create<String>("csv"));
    setValidFormats_("in_internal_calibration", ListUtils::create<String>("csv"));
	  setValidFormats_("in_protein_quantifier_protein", ListUtils::create<String>("csv"));
	  setValidFormats_("in_IDMapper", ListUtils::create<String>("FeatureXML"));
	  setValidFormats_("in_map_RT_transformer", ListUtils::create<String>("FeatureXML"));
    setValidFormats_("in_rawfiles_false_discovery_rate", ListUtils::create<String>("MzML"));
	  setValidFormats_("in_post_false_discovery_rate", ListUtils::create<String>("IdXML"));
	  setValidFormats_("in_feature_linker_unlabeledQT", ListUtils::create<String>("consensusXML"));
    setValidFormats_("in_contaminant_database", ListUtils::create<String>("Fasta"));
	  registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
  }

  ExitCodes main_(int, const char**)
  {
    StringList ins_ProteinQuantifier_Peptide = getStringList_("in_protein_quantifier_peptide");
    StringList ins_ProteinQuantifier_Protein = getStringList_("in_protein_quantifier_protein");
    StringList ins_IDMapper = getStringList_("in_IDMapper");
    StringList ins_MapRTTransformer = getStringList_("in_map_RT_transformer");
    StringList ins_rawfiles_FalseDiscoveryRate = getStringList_("in_rawfiles_false_discovery_rate");
    StringList ins_Post_FalseDiscoveryRate = getStringList_("in_post_false_discovery_rate");
    StringList ins_FeatureLinkerUnlabeledQT = getStringList_("in_feature_linker_unlabeledQT");
    StringList ins_Contaminant_DataBase = getStringList_("in_contaminant_database");
    StringList ins_Internal_Calibration = getStringList_("in_internal_calibration");
    String out = getStringOption_("out");
    QCCsvFiles qc_types_CSV;
    QCIDXMLFiles qc_types_IDXML;
    QCConsensusMaps qc_types_CONSENSUS;
    QCFeatureMaps qc_types_FEATURE;
    QCFastaFiles qc_types_FASTA;

    //Data is preprocessed for the upcoming runAllMetrics

    if (ins_ProteinQuantifier_Peptide.size()!=0)
		{
		  for(StringList::const_iterator it=ins_ProteinQuantifier_Peptide.begin();it!=ins_ProteinQuantifier_Peptide.end();++it)
			{
			  CsvFile fl(*it,'	',false,-1);
				qc_types_CSV.ProteinQuantifier_Peptide.push_back(fl);
			}
    }
    if (ins_ProteinQuantifier_Protein.size()!=0)
		{
		  for(StringList::const_iterator it=ins_ProteinQuantifier_Protein.begin();it!=ins_ProteinQuantifier_Protein.end();++it)
			{
			  CsvFile fl(*it,'	',false,-1);
				qc_types_CSV.ProteinQuantifier_Protein.push_back(fl);
			}
    }
    if (ins_MapRTTransformer.size()!=0)
		{
      vector<String> mrawfiles;
		  for(StringList::const_iterator it=ins_MapRTTransformer.begin();it!=ins_MapRTTransformer.end();++it)
			{
			  FeatureMap features;
			  FeatureXMLFile().load(*it, features);
        for(Size i = 0; i< features.size();i++)
        {
           features[i].applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
        }
        //mrawfiles.push_back(features.getMetaValue("spectra_data"));
        qc_types_FEATURE.Map_RT.push_back(features);
	    }
    }
		if (ins_IDMapper.size()!=0)
		{
      vector<String> frawfiles;
		  for(StringList::const_iterator it=ins_IDMapper.begin();it!=ins_IDMapper.end();++it)
			{
			  FeatureMap features;
			  FeatureXMLFile().load(*it, features);
        for(Size i = 0; i< features.size();i++)
        {
           features[i].applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
        }
        qc_types_FEATURE.ID_mapper.push_back(features);
	    }
    }
    if (ins_Post_FalseDiscoveryRate.size()!=0)
		{
      if(ins_rawfiles_FalseDiscoveryRate.size()!=ins_Post_FalseDiscoveryRate.size())
      {
        throw Exception::MissingInformation(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION,"invalid number of input rawfiles (rawfiles_FalseDiscoveryRate)");
      }
		  for(Size i=0;i<ins_Post_FalseDiscoveryRate.size();++i)
			{
        qc_types_IDXML.Post_False_Discovery_Rate.push_back(ins_Post_FalseDiscoveryRate[i]);
        qc_types_IDXML.Post_False_Discovery_Rate_Raw_Files.push_back(ins_rawfiles_FalseDiscoveryRate[i]);
		  }
	  }
    if(ins_FeatureLinkerUnlabeledQT.size()!=0)
		{
		  for(StringList::const_iterator it=ins_FeatureLinkerUnlabeledQT.begin();it!=ins_FeatureLinkerUnlabeledQT.end();++it)
			{
			  ConsensusMap CMap;
			  ConsensusXMLFile().load(*it,CMap);
        qc_types_CONSENSUS.Feature_Linker_Unlabled.push_back(CMap);
		  }
	  }
    if(ins_Contaminant_DataBase.size() != 0)
    {
      for(StringList::const_iterator it = ins_Contaminant_DataBase.begin() ; it != ins_Contaminant_DataBase.end(); ++it)
      {
        vector<FASTAFile::FASTAEntry> entryObj;
        FASTAFile().load(*it,entryObj);
        qc_types_FASTA.Contaminant_Database.push_back(entryObj);
      }
    }
    if(ins_Internal_Calibration.size()!=0)
    {
      for(StringList::const_iterator it = ins_Internal_Calibration.begin() ; it != ins_Internal_Calibration.end(); ++it)
      {
        CsvFile fl(*it,',',false,-1);
        qc_types_CSV.Internal_Calibration.push_back(fl);
      }
    }

    //Within class
		QCMetrics metricObj(qc_types_CSV,qc_types_FASTA,qc_types_IDXML,qc_types_FEATURE,qc_types_CONSENSUS,out);
    metricObj.runAllMetrics();
    return EXECUTION_OK;
  }
};
int main(int argc, const char** argv)
{
  QCCollector tool;
  return tool.main(argc,argv);
}
