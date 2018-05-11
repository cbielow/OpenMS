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
#include <OpenMS/CONCEPT/QCPeptideIntensity.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <stdlib.h>
#include <boost/regex.hpp>
#include <stdexcept>
#include <sstream>

using namespace OpenMS;
using namespace std;

QCPeptideIntensity::~QCPeptideIntensity()
{
}

//help function to check if the peptide is already existing in MzTab
Size QCPeptideIntensity::FindPeptideInMzTab_(const String& peptide, const vector<MzTabPeptideSectionRow>& ROWS)
{
  for(Size i = 0; i < ROWS.size(); i++)
  {
    if(peptide == ROWS[i].sequence.get())
    {
      return i;
    }
  }
  cout<<"not fond"<<endl;
  throw invalid_argument( "not found in given MzTab" );
}

bool QCPeptideIntensity::PeptideIntensity(MzTab& mztab)
{
  bool outbool = false;
  string rawfiles;
  vector<String> rawfilesdelimiter;

  //extract the peptide informations from input csvfile
  for(vector<CsvFile>::const_iterator it = CsvFilesPeptide_.begin(); it!=CsvFilesPeptide_.end();it++)
  {
   CsvFile fl = *it;
   StringList CurrentRow;
   Size line = 0;
   Size maxRow = fl.rowCount();
   MzTabPeptideSectionRows PepROWS = mztab.getPeptideSectionRows();

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
   MzTabPeptideSectionRow PepROW;
   Size peptide;
   vector<Size> abundance_list;
   Size counter= 0;
   bool headfinder = false;

   //goes through the file and extract the peptide sequences and the abundances from it.
   while(line<maxRow)
   {
     fl.getRow(line,CurrentRow);
     if(!headfinder)
     {
      for(Size j=0;j<CurrentRow.size();j++)
      {
        boost::regex rgx("abundance");
        boost::smatch match;
        if(CurrentRow[j]=="\"peptide\"")
        {
          peptide = j;
          counter++;
        }
        else if(regex_search(CurrentRow[j],match,rgx))
        {
          abundance_list.push_back(j);
          counter++;
        }
      }
      counter>1?headfinder=true:headfinder=false;
      line++;
     }
    else
    {

      //writes the peptide informations (peptide sequence and abundance) into MzTab
      try
      {
        Size mztabline = this -> FindPeptideInMzTab_(CurrentRow[peptide],PepROWS);
        MzTabPeptideSectionRow referencePepRow;
        referencePepRow = PepROWS[mztabline];
        Size pos = 0;
        for(Size v = 0; v<abundance_list.size();v++)
        {
          if(abundance_list[v]>0)
          {
            MzTabPeptideSectionRow PepROW = referencePepRow;
            MzTabString mzAbu(CurrentRow[abundance_list[v]]);
            MzTabOptionalColumnEntry abu1 = make_pair("opt_intensity",mzAbu);
            PepROW.opt_.push_back(abu1);
            outbool = true;
            if(pos == 0)
            {
              PepROWS[mztabline] = PepROW;
              pos++;
            }
            else
            {
              PepROWS.insert(PepROWS.begin()+mztabline+pos,PepROW);
              pos ++;
            }
          }
        }
      }
      //nicht sicher, ob wir den catch Fall brauchen, da PeptideIntensity die Sequenzen eigentlich nicht braucht (nach PTXQC)
      //außerdem wird PeptideIntensity immer mit ProteinAndPeptideCount ausgeführt, wodurch man die Peptide Sequenzen von dort schon erhält
      catch ( const std::invalid_argument& e )
      {
        /*
        cout<<"hier?"<<endl;
        MzTabPeptideSectionRow PepROW;
        MzTabString pep_seq(CurrentRow[peptide]);
        PepROW.sequence = pep_seq;
        MzTabString mzAbu(CurrentRow[abundance_list[0]]);
        abu1 = make_pair("opt_abundance1",mzAbu);
        PepROW.opt_.push_back(abu1);
        MzTabString mzAbu2(CurrentRow[abundance_list[0]]);
        abu2 = make_pair("opt_abundance2",mzAbu2);
        PepROW.opt_.push_back(abu2);
        MzTabString mzAbu3(CurrentRow[abundance_list[0]]);
        abu3 = make_pair("opt_abundance3", mzAbu3);
        PepROW.opt_.push_back(abu3);
        outbool = true;
        PepROWS.push_back(PepROW);
        */
      }
      line++;
    }
   }
    mztab.setPeptideSectionRows(PepROWS);
  }
return outbool;
}
