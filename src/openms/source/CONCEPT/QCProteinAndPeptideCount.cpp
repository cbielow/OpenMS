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
// $Maintainer: Mohammad El-Ismail
// $Authors: Mohammad El-Ismail
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <sstream>
#include <stdexcept>
#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

QCProteinAndPeptideCount::~QCProteinAndPeptideCount(){

}
bool QCProteinAndPeptideCount::ProtAndPepCount( MzTab& mztab)
{
MzTabPeptideSectionRows PepROWS;
MzTabProteinSectionRows ProtROWS;
string rawfiles;
vector<String> rawfilesdelimiter;

//extract the peptide informations from the input csvFile
for(vector<CsvFile>::const_iterator it = CsvFilesPeptide.begin(); it!=CsvFilesPeptide.end();it++)
{
  StringList DataList;
  bool headfinder = false;
  CsvFile fl = *it;
  Size line = 0;
  StringList CurrentRow;
  Size maxRow = fl.rowCount();
  while(fl.getRow(line,CurrentRow)==false)
  {
    line++;
  }
  //goes through the file and if the header = peptide then safe the sequences into the StringList DataList
  while(line<maxRow)
  {
    fl.getRow(line,CurrentRow);
    if(CurrentRow[0] == "\"peptide\"")
    {
      headfinder = true;
    }
    else if(headfinder==true)
    {
      DataList.push_back(CurrentRow[0]);
    }
    if(!headfinder)
    {
      line = maxRow;
    }
    line++;
  }
  //writes the peptide informations (peptide sequence & match_time_difference) into MzTab
  for(StringList::const_iterator it=DataList.begin(); it != DataList.end(); ++it)
  {
    MzTabPeptideSectionRow PepROW;
    MzTabString PepSeq;
    MzTabString MTDiff;
    MzTabOptionalColumnEntry mTOCE = make_pair("opt_Match_Time_Difference",MTDiff);
    vector<MzTabOptionalColumnEntry> optionals;
    optionals.push_back(mTOCE);
    PepROW.opt_= optionals;
    PepSeq.set(*it);
    PepROW.sequence = PepSeq;
    PepROWS.push_back(PepROW);

  }
  mztab.setPeptideSectionRows(PepROWS);
}

//extract the protein informations from the input csvFile
for(vector<CsvFile>::const_iterator it = CsvFilesProtein.begin(); it!=CsvFilesProtein.end();it++)
{
  StringList DataList;
  bool headfinder = false;
  CsvFile fl = *it;
  Size line = 0;
  StringList CurrentRow;
  Size maxRow = fl.rowCount();

  //extract the metadata rawfiles from the csvFile and puts the rawfiles into the vector rawsdelimiter
  while(fl.getRow(line,CurrentRow)==false)
  {
    boost::regex rgx("Rawfiles");
    boost::smatch match;
    bool found = boost::regex_search(CurrentRow[0],match,rgx);
    boost::regex rgx2("# Rawfiles: [[]");
    boost::regex rgx3("[]]");
    if(found)
    {
      String result2 = regex_replace(CurrentRow[0],rgx2,"");
      result2 = regex_replace(result2,rgx3,"");
      rawfiles = result2;
      if(rawfilesdelimiter.empty())
      {
        stringstream ss(rawfiles);
        while(ss.good())
        {
          String substr;
          getline(ss, substr,',');
          rawfilesdelimiter.push_back(substr);
        }
      }
    }
    line++;
  }

  //goes through the file and if the header = protein then write all the protein sequences into the StringList DataList
  while(line<maxRow)
  {
    fl.getRow(line,CurrentRow);
    if(CurrentRow[0] == "\"protein\"")
    {
      headfinder = true;
    }
    else if(headfinder==true)
    {
      DataList.push_back(CurrentRow[0]);
    }
    if(!headfinder)
    {
      line = maxRow;
    }
    line++;
  }

  //writes the protein informations (protein sequence and match_time_difference) into MzTab
  for(StringList::const_iterator it=DataList.begin(); it != DataList.end(); ++it)
  {
    MzTabProteinSectionRow ProtROW;
    MzTabString ProtSeq;
    MzTabString MTDiff;
    MzTabOptionalColumnEntry mTOCE = make_pair("opt_Match_Time_Difference",MTDiff);
    vector<MzTabOptionalColumnEntry> optionals;
    optionals.push_back(mTOCE);
    ProtROW.opt_= optionals;
    ProtSeq.set(*it);
    ProtROW.description = ProtSeq;
    ProtROWS.push_back(ProtROW);
  }
  mztab.setProteinSectionRows(ProtROWS);
}
return ProtROWS.size()==0 && PepROWS.size()== 0 ? false : true ;
}
