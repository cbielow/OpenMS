#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMS2IdentificationRate.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <boost/regex.hpp>
#include <OpenMS/FORMAT/PercolatorOutfile.h>


using namespace OpenMS;
using namespace std;

QCMS2IdentificationRate::~QCMS2IdentificationRate(){

}

int QCMS2IdentificationRate::MS2IDRateidentifier_( MzTabFile& MzTabOutputFile, String out)
{
  MzTab mztab;
  vector<pair<String,String>> idXMLFiles;
  Size control;
  boost::regex idxml("[A-Za-z0-9]+[.]idXML");
  boost::regex mzml("[A-Za-z0-9]+[.]mzML");
  boost::regex replecment("idXML");
  for(vector<pair<String,pair<String,String>>>::const_iterator it = ivec_.begin();it!=ivec_.end();++it)
  {
    if(it->first=="Post_FalseDiscoveryRate"){
      idXMLFiles.push_back(it->second);
    }
  }
  for(vector<pair<String,String>>::const_iterator it=idXMLFiles.begin();it!=idXMLFiles.end();++it){
    boost::smatch matchmzml;
    boost::smatch matchidxml;
    boost::regex_search(it->first,matchmzml,mzml);
    boost::regex_search(it->second,matchidxml,idxml);
    String rawfiles = boost::regex_replace(matchidxml[0].str(),replecment,"mzML");
    if(matchmzml[0]!=rawfiles)
    {
      throw Exception::MissingInformation(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION,"invalid order of input rawfiles_FalseDiscoveryRate or input Post_FalseDiscoveryRate Files. The Input Files must have the same order");
    }
    IdXMLFile il;
		vector<PeptideIdentification> pep_ids;
		vector<ProteinIdentification> prot_ids;
		il.load(it->second, prot_ids, pep_ids);
    MzMLFile mzmlfile;
    typedef PeakMap MapType;
    MapType exp;
    mzmlfile.getOptions().setMSLevels({2});
    mzmlfile.load(it->first,exp);
    Size MS2_spectra_count = exp.getSpectra().size();
    //cout<<"spectra MS: "<<MS2_spectra_count<<endl;
  }
}
