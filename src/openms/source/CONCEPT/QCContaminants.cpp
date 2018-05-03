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
// $Maintainer: Anton Haberland
// $Authors: Anton Haberland
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/QCContaminants.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/FASTAFile.h>

using namespace std;
using namespace OpenMS;

QCContaminants::~QCContaminants()
{

}

bool QCContaminants::QCContaminantCalculator(MzTab& mztab, bool protAndPepCount)
{
  typedef vector<pair<String,vector<FASTAFile::FASTAEntry>>> AllFastas;
  typedef vector<vector<FASTAFile::FASTAEntry>> AllEntrys;
  AllEntrys entryList;
  for(AllFastas::const_iterator it = fFiles_.begin(); it != fFiles_.end(); it ++)
  {
    if(it->first ==  "Contaminant_DataBase")
    {
       entryList.push_back(it->second);
    }
  }
  if(protAndPepCount && entryList.size() > 0)
  {
    MzTabPeptideSectionRows pepSecROWS;
    MzTabProteinSectionRows protSecROWS;
    pepSecROWS = mztab.getPeptideSectionRows();
    protSecROWS = mztab.getProteinSectionRows();
    for(MzTabPeptideSectionRows::iterator it =  pepSecROWS.begin(); it != pepSecROWS.end(); it ++)
    {
      MzTabString finder;
      finder.set("-");
      for(AllEntrys::iterator itt = entryList.begin(); itt != entryList.end(); itt ++)
      {
        for(vector<FASTAFile::FASTAEntry>::iterator entryIter = itt->begin(); entryIter != itt->end(); entryIter ++)
        {
          if(it->sequence.get() == entryIter->sequence)
          {
            finder.set("+");
            entryIter  = itt->end();
            break;
          }
        }
        if(finder.get() == "+")
        {
          break;
        }
      }
      MzTabOptionalColumnEntry col;
      col=make_pair("opt_isContaminant",finder);
      it->opt_.push_back(col);
    }
      for(MzTabProteinSectionRows::iterator it = protSecROWS.begin(); it != protSecROWS.end(); it ++)
      {
        MzTabString finder;
        finder.set("-");
        String wantedSeq = it->description.get();
        for(AllEntrys::iterator itt = entryList.begin(); itt != entryList.end(); itt ++)
        {
          for(vector<FASTAFile::FASTAEntry>::iterator entryIter = itt->begin(); entryIter != itt->end(); entryIter ++ )
          {
            if(wantedSeq == "\""+entryIter-> identifier+"\"" )
            {
              finder.set("+");
              entryIter  = itt->end();
              break;
            }
          }
          if(finder.get() == "+")
          {
            break;
          }
        }
        MzTabOptionalColumnEntry col;
        col = make_pair("opt_isContaminant",finder);
        it->opt_.push_back(col);
      }
    mztab.setPeptideSectionRows(pepSecROWS);
    mztab.setProteinSectionRows(protSecROWS);
  }
  return protAndPepCount;
}
