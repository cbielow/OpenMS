// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricKitDetection.h>
#include <OpenMS/METADATA/SILACDetector.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Detects the quantitative labelling strategy used in an LC-MS/MS experiment.

    A thin orchestrator that runs the two complementary detectors and combines their verdicts:
      - isobaric labelling (TMT / TMTpro / iTRAQ) via IsobaricKitDetection (reporter ions in MS2/MS3), and
      - metabolic SILAC labelling via SILACDetector (precursor mass-difference statistics over MS2 scans).

    These are independent — a run can be isobaric, SILAC, both, or neither (label-free) — so both verdicts
    are reported. Use detect() for the combined Result and report() for a human-readable summary.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI LabellingDetector
  {
  public:
    /// Combined detection result.
    struct OPENMS_DLLAPI Result
    {
      /// @name Isobaric (TMT / iTRAQ)
      /// @{
      bool isobaric_detected = false;                                   ///< a valid isobaric kit was detected
      IsobaricKitDetection::MethodType isobaric_kit = IsobaricKitDetection::MethodType::UNKNOWN; ///< the detected kit (UNKNOWN if none)
      std::vector<IsobaricKitDetection::KitResult> isobaric_candidates;  ///< all candidate kits, ranked by score (highest first)
      /// @}

      /// @name SILAC
      /// @{
      bool silac_applicable = false;  ///< whether SILAC detection could be run (needs MS2 precursor scans with control distances)
      bool silac_detected = false;    ///< whether a SILAC labelling was found
      SILACDetector silac;            ///< the SILAC statistics (z-scores / p-values per K4/R6-K6/K8/R10 distance)
      /// @}

      /// True if neither isobaric nor SILAC labelling was detected (i.e. likely label-free).
      bool isLabelFree() const { return !isobaric_detected && !silac_detected; }
    };

    /**
      @brief Detect isobaric and SILAC labelling in @p exp.

      Runs IsobaricKitDetection::detect() and (if the data has MS2 precursor scans) SILACDetector::detectSILAC().
      SILAC detection that cannot be run (no MS2 scans, or no control distances to normalise against) is
      reported via Result::silac_applicable = false rather than throwing.

      @param exp Input experiment.
      @param isobaric_params Thresholds forwarded to IsobaricKitDetection.
      @return The combined Result.
    */
    static Result detect(const PeakMap& exp,
                         const IsobaricKitDetection::Parameters& isobaric_params = IsobaricKitDetection::Parameters());

    /// Human-readable, multi-line summary of @p r (isobaric verdict, SILAC verdict, overall conclusion).
    static std::string report(const Result& r);
  };

} // namespace OpenMS
