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

  // Like makeExperiment(), but channel @p noisy_idx is placed with an alternating +/- @p jitter_da
  // offset each spectrum -> the channel is still matched (jitter < tolerance) but has a large Δppm scatter.
  PeakMap makeExperimentNoisy(const vector<double>& reporter_mz, Size n, Size noisy_idx, double jitter_da)
  {
    PeakMap exp;
    for (Size s = 0; s < n; ++s)
    {
      MSSpectrum spec;
      spec.setMSLevel(2);
      spec.setType(SpectrumSettings::SpectrumType::CENTROID);
      for (Size c = 0; c < reporter_mz.size(); ++c)
      {
        double mz = reporter_mz[c];
        if (c == noisy_idx) { mz += (s % 2 == 0 ? jitter_da : -jitter_da); }
        spec.emplace_back(mz, 1000.0f);
      }
      spec.emplace_back(123.5, 500.0f); // non-reporter background within the region
      spec.sortByPosition();
      exp.addSpectrum(spec);
    }
    exp.updateRanges();
    return exp;
  }

  // index of the result for a given kit
  Size kitIndex(const vector<IsobaricKitDetection::KitResult>& res, MethodType mt)
  {
    for (Size i = 0; i < res.size(); ++i) { if (res[i].type == mt) { return i; } }
    return res.size();
  }

  // find a channel by its label within a kit result (nullptr if absent)
  const IsobaricKitDetection::ChannelStats* findChannel(const IsobaricKitDetection::KitResult& kr, const std::string& name)
  {
    for (const auto& ch : kr.channels) { if (ch.name == name) { return &ch; } }
    return nullptr;
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

START_SECTION((static std::vector<ChannelRef> referenceChannels()))
{
  const auto refs = IsobaricKitDetection::referenceChannels();
  // 35 unique TMT channels (TMT35 is the TMT superset) + 8 iTRAQ channels (113..119,121), disjoint
  TEST_EQUAL(refs.size(), 43)
  // strictly ascending in m/z (and hence unique)
  bool ascending = true;
  for (Size i = 1; i < refs.size(); ++i) { if (!(refs[i].mz > refs[i - 1].mz)) { ascending = false; } }
  TEST_TRUE(ascending)
  // contains an iTRAQ-8-only channel and a TMT-35-only channel
  auto has = [&](const std::string& name) { for (const auto& r : refs) { if (r.name == name) { return true; } } return false; };
  TEST_TRUE(has("113"))    // iTRAQ 8-plex only
  TEST_TRUE(has("135CD"))  // TMT 35-plex only
}
END_SECTION

START_SECTION((static std::vector<double> channelTolerances(const std::vector<ChannelRef>& refs, double max_tolerance_ppm)))
{
  const auto refs = IsobaricKitDetection::referenceChannels();
  const auto tol = IsobaricKitDetection::channelTolerances(refs, 30.0);
  TEST_EQUAL(tol.size(), refs.size())

  double tol_114 = -1.0, mz_114 = 0.0, tol_128ND = -1.0, mz_128ND = 0.0;
  for (Size i = 0; i < refs.size(); ++i)
  {
    if (refs[i].name == "114")   { tol_114 = tol[i];   mz_114 = refs[i].mz; }    // iTRAQ: sparse neighbourhood
    if (refs[i].name == "128ND") { tol_128ND = tol[i]; mz_128ND = refs[i].mz; }  // TMTpro deuterated quartet: dense
  }
  ABORT_IF(tol_114 < 0.0 || tol_128ND < 0.0)
  // iTRAQ keeps the full 30 ppm (nearest neighbour ~1 Th away)
  TEST_REAL_SIMILAR(tol_114, 30e-6 * mz_114)
  // dense TMTpro quartet (~2.9 mDa spacing) is capped well below 30 ppm
  TEST_TRUE(tol_128ND < 30e-6 * mz_128ND)
}
END_SECTION

START_SECTION((static void classifyChannels(std::vector<ChannelStats>& channels, const std::vector<double>& channel_tol_ppm, const Parameters& params)))
{
  using ChannelStats = IsobaricKitDetection::ChannelStats;
  const IsobaricKitDetection::Parameters params; // defaults

  auto set = [](ChannelStats& c, double pop, Size n, double sd) { c.found_fraction = pop; c.n_found = n; c.ppm_spread = sd; };

  // --- presence/abundance classification + absolute noise test ---
  std::vector<ChannelStats> ch(6);
  set(ch[0], 1.0, 10, 0.5); // ok
  set(ch[1], 1.0, 10, 0.5); // ok (made noisy below)
  set(ch[2], 1.0, 10, 0.5); // ok
  set(ch[3], 1.0, 10, 0.5); // ok
  set(ch[4], 0.0,  0, 0.0); // missing
  set(ch[5], 0.05, 1, 0.0); // underpopulated (0.05 < 0.3 * 1.0)
  std::vector<double> tol_ppm(6, 12.0);

  IsobaricKitDetection::classifyChannels(ch, tol_ppm, params);
  TEST_FALSE(ch[0].is_outlier)
  TEST_TRUE(ch[4].is_outlier)
  TEST_EQUAL(ch[4].outlier_reason, "missing")
  TEST_TRUE(ch[5].is_outlier)
  TEST_EQUAL(ch[5].outlier_reason, "underpopulated")

  // absolute test: sd 6 ppm > noise_sd_frac_of_tol(0.4) * 12 ppm = 4.8 ppm
  ch[1].ppm_spread = 6.0;
  IsobaricKitDetection::classifyChannels(ch, tol_ppm, params);
  TEST_TRUE(ch[1].is_outlier)
  TEST_EQUAL(ch[1].outlier_reason, "high delta-ppm variance")
  TEST_FALSE(ch[0].is_outlier) // an exact channel stays ok

  // --- relative (median/MAD) noise test, with the absolute test disabled by a wide tolerance ---
  std::vector<ChannelStats> ch2(5);
  const double sds[] = {1.0, 2.0, 3.0, 4.0, 20.0}; // median 3, MAD 1 -> threshold 3 + 3*1 = 6
  for (Size i = 0; i < 5; ++i) { set(ch2[i], 1.0, 10, sds[i]); }
  std::vector<double> tol2(5, 200.0); // 0.4*200 = 80 ppm: absolute test cannot fire
  IsobaricKitDetection::classifyChannels(ch2, tol2, params);
  TEST_TRUE(ch2[4].is_outlier) // sd 20 > 6
  TEST_EQUAL(ch2[4].outlier_reason, "high delta-ppm variance")
  TEST_FALSE(ch2[3].is_outlier) // sd 4 < 6

  // --- median-Δppm consistency: a channel whose median offset deviates from the consensus is flagged ---
  std::vector<ChannelStats> ch3(5);
  const double offs[] = {1.0, 1.2, 0.8, 1.1, 10.0}; // consensus ~1.1; the 10 ppm channel deviates ~9 > offset_consistency_ppm (5)
  for (Size i = 0; i < 5; ++i)
  {
    set(ch3[i], 1.0, 10, 0.5);             // well-populated, low scatter (not noisy)
    ch3[i].ppm_offset = offs[i];
  }
  std::vector<double> tol3(5, 200.0);       // wide tolerance: absolute noise test cannot fire
  IsobaricKitDetection::classifyChannels(ch3, tol3, params);
  TEST_TRUE(ch3[4].is_outlier)
  TEST_EQUAL(ch3[4].outlier_reason, "median deltaPPM inconsistent")
  TEST_FALSE(ch3[0].is_outlier) // consistent channel stays ok
}
END_SECTION

START_SECTION((static double kitScore(Size num_channels, Size detected_count, Size num_covered)))
{
  TEST_REAL_SIMILAR(IsobaricKitDetection::kitScore(11, 11, 11), 1.0)         // exact match
  TEST_REAL_SIMILAR(IsobaricKitDetection::kitScore(10, 11, 10), 10.0 / 11.0) // too small (misses one detected channel)
  TEST_REAL_SIMILAR(IsobaricKitDetection::kitScore(16, 11, 11), 11.0 / 16.0) // too big (surplus channels)
  TEST_REAL_SIMILAR(IsobaricKitDetection::kitScore(4, 0, 0), 0.0)            // nothing detected
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
  TEST_EQUAL(res11.front().num_covered, 11)
  TEST_EQUAL(res11.front().num_uncovered, 0)
  // a too-small kit (TMT 10-plex) must rank below TMT 11-plex
  TEST_TRUE(res11[kitIndex(res11, MethodType::TMT_10PLEX)].score
            < res11[kitIndex(res11, MethodType::TMT_11PLEX)].score)
  // and TMT 10-plex is flagged as too small (misses 131C)
  TEST_TRUE(res11[kitIndex(res11, MethodType::TMT_10PLEX)].num_uncovered > 0)

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
  TEST_EQUAL(res6.front().num_uncovered, 0)
  // clean-signal: the 6 reporter peaks (1000 each) out of a region of 6*1000 + 500 background = 6000/6500 ~= 0.923
  TEST_TRUE(res6.front().clean_signal_fraction > 0.9 && res6.front().clean_signal_fraction <= 1.0)

  // ---- iTRAQ 4-plex sample ------------------------------------------------
  PeakMap exp_itraq = makeExperimentForKit(MethodType::ITRAQ_4PLEX, 8);
  auto res_itraq = IsobaricKitDetection::detect(exp_itraq);
  TEST_FALSE(res_itraq.empty())
  ABORT_IF(res_itraq.empty())
  TEST_TRUE(res_itraq.front().type == MethodType::ITRAQ_4PLEX)
  TEST_EQUAL(res_itraq.front().num_uncovered, 0)
  // iTRAQ 8-plex (a superset) must rank below iTRAQ 4-plex
  TEST_TRUE(res_itraq[kitIndex(res_itraq, MethodType::ITRAQ_8PLEX)].score
            < res_itraq[kitIndex(res_itraq, MethodType::ITRAQ_4PLEX)].score)

  // ---- iTRAQ keeps the full ppm tolerance (per-channel nearest-neighbour cap, not a global one) ----
  // Offset channel "115" by 2.5 mDa (~22 ppm): that exceeds the dense-TMTpro global spacing cap (~1.46 mDa)
  // but is within 30 ppm. iTRAQ channels are ~1 Th apart, so the channel must still be matched (populated).
  vector<double> itraq4_mz;
  {
    auto m = IsobaricQuantitationMethod::create(MethodType::ITRAQ_4PLEX);
    for (const auto& c : m->getChannelInformation()) { itraq4_mz.push_back(c.center); }
  }
  PeakMap exp_itraq_tol = makeExperimentNoisy(itraq4_mz, 8, 1, 0.0025);
  auto res_it_tol = IsobaricKitDetection::detect(exp_itraq_tol);
  const Size ii4 = kitIndex(res_it_tol, MethodType::ITRAQ_4PLEX);
  ABORT_IF(ii4 >= res_it_tol.size())
  const IsobaricKitDetection::ChannelStats* ch_it = findChannel(res_it_tol[ii4], "115");
  ABORT_IF(ch_it == nullptr)
  TEST_EQUAL(ch_it->n_found, 8) // would be 0 under a single global (~1.46 mDa) tolerance cap

  // ---- a well-populated but mass-inaccurate channel is flagged 'noisy' ----
  // jitter the 128C channel (index 2) by +-1.1 mDa each scan: still matched (< tolerance) but with large Δppm scatter
  PeakMap exp_noisy = makeExperimentNoisy(tmt6, 10, 2, 0.0011);
  auto res_noisy = IsobaricKitDetection::detect(exp_noisy);
  TEST_FALSE(res_noisy.empty())
  ABORT_IF(res_noisy.empty())
  const Size i6 = kitIndex(res_noisy, MethodType::TMT_6PLEX);
  ABORT_IF(i6 >= res_noisy.size())
  const IsobaricKitDetection::ChannelStats* ch_noisy = findChannel(res_noisy[i6], "128"); // TMT6 label for 128C
  const IsobaricKitDetection::ChannelStats* ch_clean = findChannel(res_noisy[i6], "126");
  ABORT_IF(ch_noisy == nullptr || ch_clean == nullptr)
  TEST_TRUE(ch_noisy->is_outlier)
  TEST_EQUAL(ch_noisy->outlier_reason, "high delta-ppm variance")
  TEST_FALSE(ch_clean->is_outlier) // an exact channel stays 'ok'

  // ---- a channel of a too-large kit that is absent from the data is flagged 'missing' ----
  const Size i11 = kitIndex(res6, MethodType::TMT_11PLEX);
  ABORT_IF(i11 >= res6.size())
  const IsobaricKitDetection::ChannelStats* ch_missing = findChannel(res6[i11], "131C"); // not part of a 6-plex sample
  ABORT_IF(ch_missing == nullptr)
  TEST_TRUE(ch_missing->is_outlier)
  TEST_EQUAL(ch_missing->outlier_reason, "missing")

  // ---- validity gate: reporter region dominated by non-reporter signal -> no valid kit (LFQ-like) ----
  // a few weak peaks at TMT positions, but a huge non-reporter peak dominates the region every scan
  PeakMap exp_lfq;
  for (Size s = 0; s < 10; ++s)
  {
    MSSpectrum spec;
    spec.setMSLevel(2);
    spec.setType(SpectrumSettings::SpectrumType::CENTROID);
    spec.emplace_back(TMTMasses::TMT_126, 100.0f);
    spec.emplace_back(TMTMasses::TMT_127N, 100.0f);
    spec.emplace_back(130.0, 100000.0f); // huge non-reporter peak (not at any channel m/z)
    spec.sortByPosition();
    exp_lfq.addSpectrum(spec);
  }
  exp_lfq.updateRanges();
  auto res_lfq = IsobaricKitDetection::detect(exp_lfq);
  TEST_FALSE(res_lfq.empty())
  ABORT_IF(res_lfq.empty())
  // no kit may be valid, and the top result must have probability 0 (nothing detected)
  TEST_FALSE(res_lfq.front().is_valid)
  TEST_REAL_SIMILAR(res_lfq.front().score, 0.0)
  for (const auto& kr : res_lfq) { TEST_FALSE(kr.is_valid) }
  // by contrast, the clean TMT 11-plex sample yields a valid top kit
  TEST_TRUE(res11.front().is_valid)

  // ---- MS3 reporters are preferred over MS2 (SPS-MS3 / MultiNotch) -------------------------------
  // MS2 carries an iTRAQ-4 pattern; MS3 carries a TMT-6 pattern. With MS3 present, only MS3 must be used.
  PeakMap exp_ms3;
  const vector<double> itraq4 = {114.1112, 115.1082, 116.1116, 117.1149};
  for (Size s = 0; s < 10; ++s)
  {
    MSSpectrum ms2;
    ms2.setMSLevel(2);
    ms2.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (double mz : itraq4) { ms2.emplace_back(mz, 1000.0f); }
    ms2.sortByPosition();
    exp_ms3.addSpectrum(ms2);

    MSSpectrum ms3;
    ms3.setMSLevel(3);
    ms3.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (double mz : tmt6) { ms3.emplace_back(mz, 1000.0f); }
    ms3.sortByPosition();
    exp_ms3.addSpectrum(ms3);
  }
  exp_ms3.updateRanges();
  auto res_ms3 = IsobaricKitDetection::detect(exp_ms3);
  TEST_FALSE(res_ms3.empty())
  ABORT_IF(res_ms3.empty())
  TEST_TRUE(res_ms3.front().type == MethodType::TMT_6PLEX)                          // MS3 (TMT6) used
  TEST_TRUE(res_ms3.front().is_valid)
  TEST_FALSE(res_ms3[kitIndex(res_ms3, MethodType::ITRAQ_4PLEX)].is_valid)          // MS2 (iTRAQ4) ignored

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
