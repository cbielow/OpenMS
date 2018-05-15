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
#pragma once
#include <vector>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <utility>
#include <regex>

namespace OpenMS
{
    struct QCFeatureMaps
    {
      std::vector<OpenMS::FeatureMap> Map_RT;
      std::vector<OpenMS::FeatureMap> ID_mapper;
      //std::vector<OpenMS::FeatureMap>
    };

    struct QCIDXMLFiles
    {
      std::vector<OpenMS::String> Post_False_Discovery_Rate_Raw_Files;
      std::vector<OpenMS::String> Post_False_Discovery_Rate;
    };

    struct QCCsvFiles
    {
      std::vector<OpenMS::CsvFile> ProteinQuantifier_Peptide;
      std::vector<OpenMS::CsvFile> ProteinQuantifier_Protein;
      std::vector<OpenMS::CsvFile> Internal_Calibration;

    };

    struct QCConsensusMaps
    {
      std::vector<OpenMS::ConsensusMap> Feature_Linker_Unlabled;
    };

    struct QCFastaFiles
    {
      std::vector<std::vector<OpenMS::FASTAFile::FASTAEntry>> Contaminant_Database;
    };


    class OPENMS_DLLAPI QCMetrics
    {


    public:
        QCMetrics(const QCCsvFiles& csv, const QCFastaFiles& fas, const QCIDXMLFiles& idx, const QCFeatureMaps& feat, const QCConsensusMaps& cons, const OpenMS::String& out):
        Feat_Maps_(feat),
        Idxml_(idx),
        CFiles_(csv),
        Consensus_Maps_(cons),
        Fasta_File_(fas),
        out_(out)
        {
        }
        ~QCMetrics();
        void runAllMetrics();
      protected:
        const QCFeatureMaps Feat_Maps_;
        const QCIDXMLFiles Idxml_;	//Peptide der IDXML's
        const QCCsvFiles CFiles_;	//Proteine der CSVFiles's;
        const QCConsensusMaps Consensus_Maps_;    //Alle ConsensusXMLFiles
        const QCFastaFiles Fasta_File_;   //Alle FastaFiles
        OpenMS::String out_;
	};
}
