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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCContaminants.h>
#include <OpenMS/CONCEPT/QCMetrics.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <OpenMS/CONCEPT/QCMS2IdentificationRate.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <vector>
#include <utility>

using namespace OpenMS;
using namespace std;

Metrics::~Metrics()
 {

 }
//Wenn ihr die Metriken schreibt lasst euch bitte ein Int ausgeben 1 wenn erfolgreich, 0  wenn nicht
//Die Daten die ihr erhaltet am besten als MetricMap, die wichtigen Funktionen dafür stehen oben
void Metrics::runAllMetrics()
{
////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
MzTabFile mzTabOutputFile;
MzTab mzTabOutput;
QCProteinAndPeptideCount ProtAndPepObj(CFiles_);
bool papc = ProtAndPepObj.ProtAndPepCount( mzTabOutput);
QCMS2IdentificationRate MS2IDRate(Idxml_);
bool mid = MS2IDRate.MS2IDRateidentifier( mzTabOutput);
QCContaminants ContaminantsObj(faFile_);
bool contam = ContaminantsObj.QCContaminantCalculator(mzTabOutput, papc);
QCMBRalignment MBRAlign(FeatMapsMBR_);
int mbra = MBRAlign.MBRAlignment( mzTabOutput);
if(papc == true){cout<<"ProteinAndPeptideCount Sucssessfull"<<endl;}
//if(mid == true){cout<<"MS2 identification Rate Sucssessfull"<<endl;}
if(contam == true){cout<<"Contaminants Sucssessfull"<<endl;}
if(mbra == 1){cout<<"MBR Alignment Sucssessfull"<<endl;}

mzTabOutputFile.store(out_,mzTabOutput);
}

////////////////Metrik2: ....................../////////////////////////////////////
////////////////Metrik3: ....................../////////////////////////////////////
////////////////Metrik4: ....................../////////////////////////////////////
////////////////Metrik5: ....................../////////////////////////////////////
