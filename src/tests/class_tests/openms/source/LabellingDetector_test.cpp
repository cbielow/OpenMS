// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/LabellingDetector.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/TMTMasses.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

namespace
{
  // TMT 11-plex experiment with @p n MS2 spectra (reporter peaks + a background peak; no precursors set)
  PeakMap tmt11Experiment(Size n)
  {
    const vector<double> tmt11 = {
      TMTMasses::TMT_126,  TMTMasses::TMT_127N, TMTMasses::TMT_127C, TMTMasses::TMT_128N,
      TMTMasses::TMT_128C, TMTMasses::TMT_129N, TMTMasses::TMT_129C, TMTMasses::TMT_130N,
      TMTMasses::TMT_130C, TMTMasses::TMT_131N, TMTMasses::TMT_131C
    };
    PeakMap exp;
    for (Size i = 0; i < n; ++i)
    {
      MSSpectrum s;
      s.setMSLevel(2);
      s.setRT(static_cast<double>(i));
      s.setType(SpectrumSettings::SpectrumType::CENTROID);
      for (double mz : tmt11) { s.emplace_back(mz, 1000.0f); }
      s.emplace_back(130.5, 500.0f);
      s.sortByPosition();
      exp.addSpectrum(s);
    }
    exp.updateRanges();
    return exp;
  }
}

START_TEST(LabellingDetector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static Result detect(const PeakMap& exp, const IsobaricKitDetection::Parameters& isobaric_params)))
{
  // ---- TMT 11-plex (isobaric detected; SILAC n/a because there are no precursors) ----
  auto r = LabellingDetector::detect(tmt11Experiment(10));
  TEST_TRUE(r.isobaric_detected)
  TEST_TRUE(r.isobaric_kit == IsobaricKitDetection::MethodType::TMT_11PLEX)
  TEST_TRUE(r.isobaric_candidates.size() > 0)
  TEST_FALSE(r.silac_applicable) // no MS2 precursor scans -> SILAC could not be run (handled gracefully, no throw)
  TEST_FALSE(r.silac_detected)
  TEST_FALSE(r.isLabelFree())

  // ---- label-free (a single MS1 spectrum) ----
  PeakMap lf;
  MSSpectrum ms1;
  ms1.setMSLevel(1);
  ms1.emplace_back(500.0, 1000.0f);
  lf.addSpectrum(ms1);
  auto r_lf = LabellingDetector::detect(lf);
  TEST_FALSE(r_lf.isobaric_detected)
  TEST_FALSE(r_lf.silac_detected)
  TEST_TRUE(r_lf.isLabelFree())
}
END_SECTION

START_SECTION((static Result detect(const FeatureMap& features)))
{
  FeatureMap fm;
  for (Size i = 0; i < 5; ++i)
  {
    Feature f;
    f.setRT(static_cast<double>(i));
    f.setMZ(500.0 + i);
    f.setCharge(2);
    fm.push_back(f);
  }
  auto r = LabellingDetector::detect(fm);
  TEST_FALSE(r.isobaric_applicable) // featureXML has no MS2 reporter-ion spectra
  TEST_FALSE(r.isobaric_detected)
  TEST_TRUE(r.silac_applicable)     // SILAC could be run (features present)

  // empty map -> SILAC not applicable, isobaric not applicable
  auto r_empty = LabellingDetector::detect(FeatureMap());
  TEST_FALSE(r_empty.isobaric_applicable)
  TEST_FALSE(r_empty.silac_applicable)
  TEST_TRUE(r_empty.isLabelFree())
}
END_SECTION

START_SECTION((static Result detect(const ConsensusMap& consensus)))
{
  ConsensusMap cm;
  for (Size i = 0; i < 5; ++i)
  {
    ConsensusFeature cf;
    cf.setRT(static_cast<double>(i));
    cf.setMZ(500.0 + i);
    cf.setCharge(2);
    Feature sub;
    sub.setRT(static_cast<double>(i));
    sub.setMZ(500.0 + i);
    sub.setCharge(2);
    cf.insert(0, sub); // one subfeature -> one data point
    cm.push_back(cf);
  }
  auto r = LabellingDetector::detect(cm);
  TEST_FALSE(r.isobaric_applicable) // consensusXML has no MS2 reporter-ion spectra
  TEST_FALSE(r.isobaric_detected)
  TEST_TRUE(r.silac_applicable)     // SILAC could be run (consensus features present)

  // empty map -> nothing applicable
  auto r_empty = LabellingDetector::detect(ConsensusMap());
  TEST_FALSE(r_empty.isobaric_applicable)
  TEST_FALSE(r_empty.silac_applicable)
  TEST_TRUE(r_empty.isLabelFree())
}
END_SECTION

START_SECTION((static std::string report(const Result& r)))
{
  auto r = LabellingDetector::detect(tmt11Experiment(10));
  const std::string rep = LabellingDetector::report(r);
  TEST_EQUAL(rep.find("Labelling detection") != std::string::npos, true)
  TEST_EQUAL(rep.find("TMT 11-plex") != std::string::npos, true)
  TEST_EQUAL(rep.find("SILAC") != std::string::npos, true)

  // feature/consensus input: isobaric is reported as not applicable
  FeatureMap fm;
  Feature f;
  f.setRT(1.0);
  f.setMZ(500.0);
  f.setCharge(2);
  fm.push_back(f);
  const std::string rep_feat = LabellingDetector::report(LabellingDetector::detect(fm));
  TEST_EQUAL(rep_feat.find("n/a (requires MS2/MS3 reporter-ion spectra)") != std::string::npos, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
