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
// $Maintainer: Anton Haberland, Leo Wurth, Mohammad El-Ismail$
// $Authors: Anton Haberland, Leo Wurth, Mohammad El-Ismail $
// --------------------------------------------------------------------------
#include <vector>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>
#include <regex>


    class OPENMS_DLLAPI Metrics
    {
    //das sind die eingelesenen Daten, sie können von den Metriken gelesen aber nicht umgeschrieben werden//
    //  die Vectoren enthalten alle eingelesenen Datein in den entsperechenden Formaten. Der erste Wert in jedem Paar
    //  gibt an aus welchem TOPPTOOL die Datei kommt.

    public:
        Metrics(const std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> fvec, const std::vector<std::pair<OpenMS::String,OpenMS::String>> ivec, const std::vector<std::pair<OpenMS::String,OpenMS::CsvFile>> cvecs,const std::vector<std::pair<OpenMS::String,OpenMS::ConsensusMap>> CMapVec, OpenMS::String out):
        FeatMaps_(fvec),
        Idxml_(ivec),
        CFiles_(cvecs),
        ConsensusMaps_(CMapVec),
        out_(out)
        {
        }
        ~Metrics();
        void runAllMetrics();
      protected:
        //int ProteinAndPeptideCount_(MetricMap& outPep,MetricMap& outProt) const;
        const std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> FeatMaps_;			//Alle FeatureXML Datein
        const std::vector<std::pair<OpenMS::String,OpenMS::String>> Idxml_;	//Peptide der IDXML's
        const std::vector<std::pair<OpenMS::String,OpenMS::CsvFile>> CFiles_;	//alle CSVFILES;
        const std::vector<std::pair<OpenMS::String,OpenMS::ConsensusMap>> ConsensusMaps_;    //Alle ConsensusXMLFiles
        OpenMS::String out_;
	//->Hier Metriken deklarieren<-//
	};



    	//->Hier die Metriken definieren<-//
