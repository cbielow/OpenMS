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
#include <OpenMS/CONCEPT/QCMSRecalibrationerror.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <boost/regex.hpp>
#include <vector>


using namespace OpenMS;
using namespace std;

QCMSRecalibrationerror::~QCMSRecalibrationerror(){

}

bool QCMSRecalibrationerror::QCMSRecalerror( MzTab& mztab)
{
    bool outbool = false;
    for(vector<CsvFile>::const_iterator it=cvec_.begin(); it!=cvec_.end(); ++it)
    {
      CsvFile fl = *it;
      Size maxRow = fl.rowCount();
      Size line = 0;
      StringList CurrentRow;
      while(fl.getRow(line,CurrentRow)==false)
      {
        line++;
      }
      Size intensity;
      Size mzref;
      Size mzbefore;
      Size mzafter;
      Size counter = 0;
      Internal_Calibration_Data internal_cal_data;
      bool headfinder = false;
      while(line<maxRow)
      {
        fl.getRow(line,CurrentRow);
        if(!headfinder)
        {
          for(Size i=0;i<CurrentRow.size();i++)
          {
            if(CurrentRow[i]==" intensity")
            {
              intensity = i;
              counter++;
            }
            else if(CurrentRow[i]==" mz ref")
            {
              mzref = i;
              counter++;
            }
            else if(CurrentRow[i]==" mz before")
            {
              mzbefore = i;
              counter++;
            }
            else if(CurrentRow[i]==" mz after")
            {
              mzafter = i;
              counter++;
            }
          }
          counter==4?headfinder=true:headfinder=false;
          line++;
        }
        else
        {
          internal_cal_data.intensityvec.push_back(CurrentRow[intensity]);
          internal_cal_data.mzrefvec.push_back(CurrentRow[mzref]);
          internal_cal_data.mzbeforevec.push_back(CurrentRow[mzbefore]);
          internal_cal_data.mzaftervec.push_back(CurrentRow[mzafter]);
          line++;
        }
      }
      IDMapper idm;
      FeatureMap feat;
      vector<PeptideIdentification> pepID;
      vector<ProteinIdentification> protID;
      idm.annotate(feat,pepID,protID,true,false);
    }
    return outbool;
}
