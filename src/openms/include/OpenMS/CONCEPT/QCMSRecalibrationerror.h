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
#ifndef OPENMS_CONCEPT_QCMSRecalibrationerror_H
#define OPENMS_CONCEPT_QCMSRecalibrationerror_H

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <vector>
#include <utility>

class OPENMS_DLLAPI QCMSRecalibrationerror
{
std::vector<OpenMS::CsvFile> cvec_;
  public:
    QCMSRecalibrationerror(std::vector<OpenMS::CsvFile> files):
      cvec_(files)
      {
      }
      ~QCMSRecalibrationerror();
      bool QCMSRecalerror(OpenMS::MzTab& mztab);
};
struct OPENMS_DLLAPI Internal_Calibration_Data
{
  std::vector<OpenMS::String> intensityvec;
  std::vector<OpenMS::String> mzrefvec;
  std::vector<OpenMS::String> mzbeforevec;
  std::vector<OpenMS::String> mzaftervec;
};
#endif
