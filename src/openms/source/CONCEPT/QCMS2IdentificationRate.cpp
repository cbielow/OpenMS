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
// $Maintainer: Anton Haberland, Leo Wurth$
// $Authors: Anton Haberland, Leo Wurth$
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMS2IdentificationRate.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <boost/regex.hpp>
#include <OpenMS/FORMAT/PercolatorOutfile.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <vector>


using namespace OpenMS;
using namespace std;

QCMS2IdentificationRate::~QCMS2IdentificationRate(){

}

bool QCMS2IdentificationRate::MS2IDRateidentifier( MzTab& mztab)
{
  bool outbool = false;
  vector<pair<String,String>> idXMLFiles;
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
    MapType dummy;
    MapType exp;
    mzmlfile.getOptions().setMSLevels({2});
    Size scount;
    Size ccount;
    mzmlfile.loadSize(it->first,scount,ccount);
    double identification_rate = (double)pep_ids.size()/scount;
    //////////////////////////////////////////////////////////
    MzTabParameter MetaRawAndId;
    String theRate(to_string(identification_rate));
    MetaRawAndId.setCVLabel("RawFile:");
    MetaRawAndId.setAccession(rawfiles);
    MetaRawAndId.setName("IdentificationRate:");
    MetaRawAndId.setValue(theRate);
    MzTabMetaData mzmeta = mztab.getMetaData();
    map<Size,MzTabParameter> savedMetaData = mzmeta.custom;
    savedMetaData[savedMetaData.size()+1]= MetaRawAndId;
    mzmeta.custom = savedMetaData;
    mztab.setMetaData(mzmeta);
    outbool = true;

  }
  return outbool;
}
