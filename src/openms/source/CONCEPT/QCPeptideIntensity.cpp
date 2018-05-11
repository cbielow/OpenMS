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

using namespace OpenMS;
using namespace std;

QCPeptideIntensity::~QCPeptideIntensity()
{
}

Size QCPeptideIntensity::FindPeptideInMzTab_(const String& peptide, const MzTab& mztab)
{
  MzTabPeptideSectionRows ROWS = mztab.getPeptideSectionRows();
  for(Size i = 0; i < ROWS.size(); i++)
  {
    if(peptide == ROWS[i].sequence.get())
    {
      return i;
    }
  }
  throw invalid_argument( "not found in given MzTab" );
}

bool QCPeptideIntensity::PeptideIntensity(MzTab& mztab)
{
  bool outbool = false;
  for(vector<CsvFile>::const_iterator it = CsvFilesPeptide_.begin(); it!=CsvFilesPeptide_.end();it++)
  {
   CsvFile fl = *it;
   StringList CurrentRow;
   Size line = 0;
   Size maxRow = fl.rowCount();
   MzTabPeptideSectionRows PepROWS = mztab.getPeptideSectionRows();

   while(fl.getRow(line,CurrentRow)==false)
   {
     line++;
   }
   MzTabPeptideSectionRow PepROW;
   MzTabOptionalColumnEntry abu1;
   MzTabOptionalColumnEntry abu2;
   MzTabOptionalColumnEntry abu3;
   Size peptide;
   Size abundance_1;
   Size abundance_2;
   Size abundance_3;
   Size counter= 0;
   bool headfinder = false;
   while(line<maxRow)
   {
     fl.getRow(line,CurrentRow);
     if(!headfinder)
     {
      for(Size i=0;i<CurrentRow.size();i++)
      {
        if(CurrentRow[i]=="\"peptide\"")
        {
          peptide = i;
          counter++;
        }
        if(CurrentRow[i]=="\"abundance_1\"")
        {
          abundance_1 = i;
          counter++;
        }
        if(CurrentRow[i]=="\"abundance_2\"")
        {
          abundance_2 = i;
          counter++;
        }
        if(CurrentRow[i]=="\"abundance_3\"")
        {
          abundance_3 = i;
          counter++;
        }
      }
      counter==4?headfinder=true:headfinder=false;
      line++;
     }

    try
    {
      Size mztabline = this -> FindPeptideInMzTab_(CurrentRow[peptide],mztab);
      MzTabPeptideSectionRow PepROW;
      PepROW = PepROWS[mztabline];
      MzTabString mzAbu(CurrentRow[abundance_1]);
      abu1 = make_pair("opt_abundance1",mzAbu);
      PepROW.opt_.push_back(abu1);
      MzTabString mzAbu2(CurrentRow[abundance_2]);
      abu2 = make_pair("opt_abundance2",mzAbu2);
      PepROW.opt_.push_back(abu2);
      MzTabString mzAbu3(CurrentRow[abundance_3]);
      abu3 = make_pair("opt_abundance3", mzAbu3);
      PepROW.opt_.push_back(abu3);
      outbool = true;
      PepROWS[mztabline] = PepROW;
    }
    catch ( const std::invalid_argument& e )
    {
      MzTabPeptideSectionRow PepROW;
      MzTabString pep_seq(CurrentRow[peptide]);
      PepROW.sequence = pep_seq;
      MzTabString mzAbu(CurrentRow[abundance_1]);
      abu1 = make_pair("opt_abundance1",mzAbu);
      PepROW.opt_.push_back(abu1);
      MzTabString mzAbu2(CurrentRow[abundance_2]);
      abu2 = make_pair("opt_abundance2",mzAbu2);
      PepROW.opt_.push_back(abu2);
      MzTabString mzAbu3(CurrentRow[abundance_3]);
      abu3 = make_pair("opt_abundance3", mzAbu3);
      PepROW.opt_.push_back(abu3);
      outbool = true;
      PepROWS.push_back(PepROW);
    }
    line++;
   }
  mztab.setPeptideSectionRows(PepROWS);
  }
return outbool;
}
