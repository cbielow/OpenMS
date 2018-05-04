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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetrics.h>
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
	  registerInputFileList_("in_ProteinQuantifierPeptide","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_ProteinQuantifierProtein","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_IDMapper","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_MapRTTransformer","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_rawfiles_FalseDiscoveryRate","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_Post_FalseDiscoveryRate","<files>", StringList(), "Input files",false,false);
	  registerInputFileList_("in_FeatureLinkerUnlabeledQT","<files>", StringList(), "Input files",false,false);
    registerInputFileList_("in_Contaminant_DataBase","<files>",StringList(), "Input files",false,false);
    setValidFormats_("in_ProteinQuantifierPeptide", ListUtils::create<String>("csv"));
	  setValidFormats_("in_ProteinQuantifierProtein", ListUtils::create<String>("csv"));
	  setValidFormats_("in_IDMapper", ListUtils::create<String>("FeatureXML"));
	  setValidFormats_("in_MapRTTransformer", ListUtils::create<String>("FeatureXML"));
    setValidFormats_("in_rawfiles_FalseDiscoveryRate", ListUtils::create<String>("MzML"));
	  setValidFormats_("in_Post_FalseDiscoveryRate", ListUtils::create<String>("IdXML"));
	  setValidFormats_("in_FeatureLinkerUnlabeledQT", ListUtils::create<String>("consensusXML"));
    setValidFormats_("in_Contaminant_DataBase", ListUtils::create<String>("Fasta"));
	  registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
  }

  ExitCodes main_(int, const char**)
  {
    StringList ins_ProteinQuantifier_Peptide = getStringList_("in_ProteinQuantifierPeptide");
    StringList ins_ProteinQuantifier_Protein = getStringList_("in_ProteinQuantifierProtein");
    StringList ins_IDMapper = getStringList_("in_IDMapper");
    StringList ins_MapRTTransformer = getStringList_("in_MapRTTransformer");
    StringList ins_rawfiles_FalseDiscoveryRate = getStringList_("in_rawfiles_FalseDiscoveryRate");
    StringList ins_Post_FalseDiscoveryRate = getStringList_("in_Post_FalseDiscoveryRate");
    StringList ins_FeatureLinkerUnlabeledQT = getStringList_("in_FeatureLinkerUnlabeledQT");
    StringList ins_Contaminant_DataBase = getStringList_("in_Contaminant_DataBase");
    String out = getStringOption_("out");
    vector<pair<String,FeatureMap>> fvec;
    vector<pair<String,FeatureMap>> mbravec;
    vector<pair<String,CsvFile>> cvec;
    vector <pair<String,ConsensusMap>> CMapVec;
    vector<pair<String,pair<String,String>>> ivec;
    vector<pair<String,vector<FASTAFile::FASTAEntry>>> faFiles;
    if (ins_ProteinQuantifier_Peptide.size()!=0)
		{
		  for(StringList::const_iterator it=ins_ProteinQuantifier_Peptide.begin();it!=ins_ProteinQuantifier_Peptide.end();++it)
			{
			  CsvFile fl(*it,'	',false,-1);
				cvec.push_back(make_pair("ProteinQuantifier_Peptide",fl));
			}
    }
    if (ins_ProteinQuantifier_Protein.size()!=0)
		{
		  for(StringList::const_iterator it=ins_ProteinQuantifier_Protein.begin();it!=ins_ProteinQuantifier_Protein.end();++it)
			{
			  CsvFile fl(*it,'	',false,-1);
				cvec.push_back(make_pair("ProteinQuantifier_Protein",fl));
			}
    }
    if (ins_MapRTTransformer.size()!=0)
		{
      vector<String> mrawfiles;
		  for(StringList::const_iterator it=ins_MapRTTransformer.begin();it!=ins_MapRTTransformer.end();++it)
			{
			  FeatureMap features;
			  FeatureXMLFile().load(*it, features);
        mrawfiles.push_back(features.getMetaValue("spectra_data"));
        mbravec.push_back(make_pair("MapRTTransformer",features));
	    }
    }
		if (ins_IDMapper.size()!=0)
		{
      vector<String> frawfiles;
		  for(StringList::const_iterator it=ins_IDMapper.begin();it!=ins_IDMapper.end();++it)
			{
			  FeatureMap features;
			  FeatureXMLFile().load(*it, features);
        frawfiles.push_back(features.getMetaValue("spectra_data"));
        fvec.push_back(make_pair("IDMapper",features));
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
        cout<<"ins Post"<<ins_Post_FalseDiscoveryRate[i]<<endl;
        cout<<"ins rawfiles"<<ins_rawfiles_FalseDiscoveryRate[i]<<endl;
			  ivec.push_back(make_pair("Post_FalseDiscoveryRate",make_pair(ins_rawfiles_FalseDiscoveryRate[i],ins_Post_FalseDiscoveryRate[i])));
		  }
	  }
    if(ins_FeatureLinkerUnlabeledQT.size()!=0)
		{
      vector<String> crawfiles;
		  for(StringList::const_iterator it=ins_FeatureLinkerUnlabeledQT.begin();it!=ins_FeatureLinkerUnlabeledQT.end();++it)
			{
			  ConsensusMap CMap;
			  ConsensusXMLFile().load(*it,CMap);
        crawfiles.push_back(CMap.getMetaValue("spectra_data"));
			  CMapVec.push_back(make_pair("FeatureLinkerUnlabeledQT",CMap));
		  }
	  }
    if(ins_Contaminant_DataBase.size() != 0)
    {
      for(StringList::const_iterator it = ins_Contaminant_DataBase.begin() ; it != ins_Contaminant_DataBase.end(); ++it)
      {
        vector<FASTAFile::FASTAEntry> entryObj;
        FASTAFile().load(*it,entryObj);
        faFiles.push_back( make_pair("Contaminant_DataBase",entryObj) );
      }
    }
		Metrics metricObj(fvec,mbravec,ivec,cvec,CMapVec,faFiles,out);

	    metricObj.runAllMetrics();
  return EXECUTION_OK;
  }
};
int main(int argc, const char** argv)
{
  QCCollector tool;
  return tool.main(argc,argv);
}
