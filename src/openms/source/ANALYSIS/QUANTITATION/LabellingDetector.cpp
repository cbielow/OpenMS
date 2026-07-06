// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/LabellingDetector.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <sstream>

namespace OpenMS
{
  LabellingDetector::Result LabellingDetector::detect(const PeakMap& exp, const IsobaricKitDetection::Parameters& isobaric_params)
  {
    Result r;

    // 1) isobaric labelling (TMT / TMTpro / iTRAQ): reporter ions in MS2/MS3
    r.isobaric_candidates = IsobaricKitDetection::detect(exp, isobaric_params);
    if (!r.isobaric_candidates.empty() && r.isobaric_candidates.front().is_valid && r.isobaric_candidates.front().score > 0.0)
    {
      r.isobaric_detected = true;
      r.isobaric_kit = r.isobaric_candidates.front().type;
    }

    // 2) SILAC: precursor mass-difference statistics over the MS2 scans.
    // SILACDetector no longer throws (empty/unsorted/no-control are handled internally); a SILAC verdict
    // can be computed whenever there is usable MS2 precursor data (MS2 scans that carry a precursor).
    if (!exp.empty() && exp.containsScanOfLevel(2))
    {
      const auto ms2 = r.silac.msExperimentToMS2Data(exp);
      if (!ms2.empty())
      {
        r.silac_applicable = true;
        r.silac_detected = r.silac.detectSILAC(ms2);
      }
    }

    return r;
  }

  LabellingDetector::Result LabellingDetector::detect(const FeatureMap& features)
  {
    Result r;
    // isobaric detection needs MS2/MS3 reporter-ion spectra, which a FeatureMap does not carry
    r.isobaric_applicable = false;

    // SILAC: precursor mass-difference statistics, one data point per feature
    const auto ms2 = r.silac.featureMapToMS2Data(features);
    if (!ms2.empty())
    {
      r.silac_applicable = true;
      r.silac_detected = r.silac.detectSILAC(ms2);
    }
    return r;
  }

  LabellingDetector::Result LabellingDetector::detect(const ConsensusMap& consensus)
  {
    Result r;
    // isobaric detection needs MS2/MS3 reporter-ion spectra, which a ConsensusMap does not carry
    r.isobaric_applicable = false;

    // SILAC: precursor mass-difference statistics, one data point per subfeature (or consensus feature if it has none)
    const auto ms2 = r.silac.consensusMapToMS2Data(consensus);
    if (!ms2.empty())
    {
      r.silac_applicable = true;
      r.silac_detected = r.silac.detectSILAC(ms2);
    }
    return r;
  }

  std::string LabellingDetector::report(const Result& r)
  {
    std::ostringstream os;
    os << "-- Labelling detection --\n";

    os << "Isobaric (TMT/iTRAQ): ";
    if (!r.isobaric_applicable) { os << "n/a (requires MS2/MS3 reporter-ion spectra)"; }
    else if (r.isobaric_detected) { os << IsobaricKitDetection::methodName(r.isobaric_kit); }
    else { os << "not detected"; }
    os << '\n';

    os << "SILAC: ";
    if (!r.silac_applicable) { os << "n/a (no MS2 precursor scans / no control distances)"; }
    else if (r.silac_detected) { os << "detected"; }
    else { os << "not detected"; }
    os << '\n';
    if (r.silac_applicable) { os << r.silac; } // detailed z-scores / p-values per distance

    os << "\nConclusion: ";
    if (r.isobaric_detected && r.silac_detected)
    {
      os << "combined isobaric (" << IsobaricKitDetection::methodName(r.isobaric_kit) << ") + SILAC labelling";
    }
    else if (r.isobaric_detected) { os << IsobaricKitDetection::methodName(r.isobaric_kit); }
    else if (r.silac_detected) { os << "SILAC"; }
    else { os << "no labelling detected (likely label-free)"; }
    os << '\n';

    return os.str();
  }

} // namespace OpenMS
