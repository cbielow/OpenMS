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
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricKitDetection.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/TMTMasses.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

using MethodType = IsobaricQuantitationMethod::MethodType;

namespace
{
  // Build a synthetic experiment with @p n MS2 spectra, each carrying the given reporter ions
  // (plus one non-reporter background peak in the reporter region).
  PeakMap makeExperiment(const vector<double>& reporter_mz, Size n)
  {
    PeakMap exp;
    for (Size s = 0; s < n; ++s)
    {
      MSSpectrum spec;
      spec.setMSLevel(2);
      spec.setType(SpectrumSettings::SpectrumType::CENTROID);
      for (double mz : reporter_mz)
      {
        spec.emplace_back(mz, 1000.0f);
      }
      spec.emplace_back(123.5, 500.0f); // non-reporter background within the region
      spec.sortByPosition();
      exp.addSpectrum(spec);
    }
    exp.updateRanges();
    return exp;
  }

  // Build an experiment for a given kit using its factory-defined channel m/z.
  PeakMap makeExperimentForKit(MethodType mt, Size n)
  {
    vector<double> mzs;
    auto method = IsobaricQuantitationMethod::create(mt);
    for (const auto& c : method->getChannelInformation()) { mzs.push_back(c.center); }
    return makeExperiment(mzs, n);
  }

  // index of the result for a given kit
  Size kitIndex(const vector<IsobaricKitDetection::KitResult>& res, MethodType mt)
  {
    for (Size i = 0; i < res.size(); ++i) { if (res[i].type == mt) { return i; } }
    return res.size();
  }
}

START_TEST(IsobaricKitDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static const std::vector<MethodType>& supportedKits()))
{
  const auto& kits = IsobaricKitDetection::supportedKits();
  // all concrete IsobaricQuantitationMethod::MethodType values (TMT 6/10/11/16/18/32/35 + iTRAQ 4/8)
  TEST_EQUAL(kits.size(), 9)
  TEST_TRUE(std::find(kits.begin(), kits.end(), MethodType::TMT_6PLEX) != kits.end())
  TEST_TRUE(std::find(kits.begin(), kits.end(), MethodType::TMT_35PLEX) != kits.end())
  TEST_TRUE(std::find(kits.begin(), kits.end(), MethodType::ITRAQ_4PLEX) != kits.end())
  TEST_TRUE(std::find(kits.begin(), kits.end(), MethodType::ITRAQ_8PLEX) != kits.end())
  // UNKNOWN must not be considered a kit
  TEST_FALSE(std::find(kits.begin(), kits.end(), MethodType::UNKNOWN) != kits.end())
}
END_SECTION

START_SECTION((static std::string methodName(MethodType mt)))
{
  TEST_EQUAL(IsobaricKitDetection::methodName(MethodType::TMT_16PLEX),
             std::string(IsobaricQuantitationMethod::methodDisplayName(MethodType::TMT_16PLEX)))
  TEST_EQUAL(IsobaricKitDetection::methodName(MethodType::ITRAQ_4PLEX),
             std::string(IsobaricQuantitationMethod::methodDisplayName(MethodType::ITRAQ_4PLEX)))
}
END_SECTION

START_SECTION((static std::vector<HierarchyNode> buildHierarchy()))
{
  const auto nodes = IsobaricKitDetection::buildHierarchy();
  TEST_EQUAL(nodes.size(), 9)

  auto find = [&](MethodType mt) -> const IsobaricKitDetection::HierarchyNode&
  {
    for (const auto& n : nodes) { if (n.type == mt) { return n; } }
    throw "not found";
  };

  // TMT 6-plex is a subset of TMT 10-plex (and has no smaller subset)
  const auto& n6 = find(MethodType::TMT_6PLEX);
  TEST_EQUAL(n6.children.size(), 0)
  TEST_TRUE(std::find(n6.parents.begin(), n6.parents.end(), MethodType::TMT_10PLEX) != n6.parents.end())

  // TMT 35-plex is the global TMT superset -> no parents
  const auto& n35 = find(MethodType::TMT_35PLEX);
  TEST_EQUAL(n35.parents.size(), 0)

  // TMT 16-plex is a subset of both 18-plex and 32-plex
  const auto& n16 = find(MethodType::TMT_16PLEX);
  TEST_TRUE(std::find(n16.parents.begin(), n16.parents.end(), MethodType::TMT_18PLEX) != n16.parents.end())
  TEST_TRUE(std::find(n16.parents.begin(), n16.parents.end(), MethodType::TMT_32PLEX) != n16.parents.end())

  // iTRAQ 4-plex is a subset of iTRAQ 8-plex (and disjoint from the TMT tree)
  const auto& ni4 = find(MethodType::ITRAQ_4PLEX);
  TEST_EQUAL(ni4.children.size(), 0)
  TEST_TRUE(std::find(ni4.parents.begin(), ni4.parents.end(), MethodType::ITRAQ_8PLEX) != ni4.parents.end())
  TEST_FALSE(std::find(ni4.parents.begin(), ni4.parents.end(), MethodType::TMT_6PLEX) != ni4.parents.end())

  const auto& ni8 = find(MethodType::ITRAQ_8PLEX);
  TEST_EQUAL(ni8.parents.size(), 0)
  TEST_TRUE(std::find(ni8.children.begin(), ni8.children.end(), MethodType::ITRAQ_4PLEX) != ni8.children.end())
}
END_SECTION

START_SECTION((static std::vector<KitResult> detect(const PeakMap& exp, const Parameters& params)))
{
  // ---- TMT 11-plex sample -------------------------------------------------
  const vector<double> tmt11 = {
    TMTMasses::TMT_126,  TMTMasses::TMT_127N, TMTMasses::TMT_127C, TMTMasses::TMT_128N,
    TMTMasses::TMT_128C, TMTMasses::TMT_129N, TMTMasses::TMT_129C, TMTMasses::TMT_130N,
    TMTMasses::TMT_130C, TMTMasses::TMT_131N, TMTMasses::TMT_131C
  };
  PeakMap exp11 = makeExperiment(tmt11, 10);
  auto res11 = IsobaricKitDetection::detect(exp11);
  TEST_FALSE(res11.empty())
  ABORT_IF(res11.empty())
  // most likely kit is TMT 11-plex and it explains all present channels
  TEST_TRUE(res11.front().type == MethodType::TMT_11PLEX)
  TEST_EQUAL(res11.front().num_explained, 11)
  TEST_EQUAL(res11.front().num_unexplained_present, 0)
  // a too-small kit (TMT 10-plex) must rank below TMT 11-plex
  TEST_TRUE(res11[kitIndex(res11, MethodType::TMT_10PLEX)].probability
            < res11[kitIndex(res11, MethodType::TMT_11PLEX)].probability)
  // and TMT 10-plex is flagged as too small (misses 131C)
  TEST_TRUE(res11[kitIndex(res11, MethodType::TMT_10PLEX)].num_unexplained_present > 0)

  // ---- TMT 6-plex sample --------------------------------------------------
  const vector<double> tmt6 = {
    TMTMasses::TMT_126, TMTMasses::TMT_127N, TMTMasses::TMT_128C,
    TMTMasses::TMT_129N, TMTMasses::TMT_130C, TMTMasses::TMT_131N
  };
  PeakMap exp6 = makeExperiment(tmt6, 8);
  auto res6 = IsobaricKitDetection::detect(exp6);
  TEST_FALSE(res6.empty())
  ABORT_IF(res6.empty())
  TEST_TRUE(res6.front().type == MethodType::TMT_6PLEX)
  TEST_EQUAL(res6.front().num_unexplained_present, 0)

  // ---- iTRAQ 4-plex sample ------------------------------------------------
  PeakMap exp_itraq = makeExperimentForKit(MethodType::ITRAQ_4PLEX, 8);
  auto res_itraq = IsobaricKitDetection::detect(exp_itraq);
  TEST_FALSE(res_itraq.empty())
  ABORT_IF(res_itraq.empty())
  TEST_TRUE(res_itraq.front().type == MethodType::ITRAQ_4PLEX)
  TEST_EQUAL(res_itraq.front().num_unexplained_present, 0)
  // iTRAQ 8-plex (a superset) must rank below iTRAQ 4-plex
  TEST_TRUE(res_itraq[kitIndex(res_itraq, MethodType::ITRAQ_8PLEX)].probability
            < res_itraq[kitIndex(res_itraq, MethodType::ITRAQ_4PLEX)].probability)

  // ---- no reporter signal -> empty result --------------------------------
  PeakMap empty_exp;
  MSSpectrum ms1;
  ms1.setMSLevel(1);
  ms1.emplace_back(500.0, 1000.0f);
  empty_exp.addSpectrum(ms1);
  auto res_none = IsobaricKitDetection::detect(empty_exp);
  TEST_TRUE(res_none.empty())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
