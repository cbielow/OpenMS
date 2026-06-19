// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricKitDetection.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <set>

namespace OpenMS
{
  namespace
  {
    /// Stable integer key for a reporter m/z, used to relate identical channels across kits.
    /// Exact bit pattern of the double: the same physical channel is populated from the same shared
    /// constant in every kit (see TMTMasses.h / the iTRAQ method literals), so identical channels yield
    /// bit-identical doubles and therefore equal keys. For positive, finite doubles the bit pattern is
    /// monotonic in value, so a std::map keyed on it also iterates in ascending m/z order.
    inline uint64_t massKey(double mz)
    {
      static_assert(sizeof(uint64_t) == sizeof(double), "massKey assumes double and uint64_t have the same size");
      uint64_t dest;
      std::memcpy(&dest, &mz, sizeof dest);
      return dest;
    }

    /// Channel labels and centers of one kit (created via the IsobaricQuantitationMethod factory).
    std::vector<std::pair<std::string, double>> getKitChannels(IsobaricQuantitationMethod::MethodType mt)
    {
      std::vector<std::pair<std::string, double>> out;
      auto method = IsobaricQuantitationMethod::create(mt);
      if (!method) { return out; }
      for (const auto& c : method->getChannelInformation())
      {
        out.emplace_back(c.name, c.center);
      }
      return out;
    }

    /// median of a (small) value vector; returns 0 for an empty vector. May reorder @p v.
    double medianOf(std::vector<double>& v)
    {
      if (v.empty()) { return 0.0; }
      return Math::median(v.begin(), v.end(), false);
    }

    /// Per-kit layout for the validity gate: the kit's reporter region and the indices of its channels in @p refs.
    struct KitLayout
    {
      double region_lo = 0.0;        ///< kit's lowest channel m/z minus the buffer
      double region_hi = 0.0;        ///< kit's highest channel m/z plus the buffer
      std::vector<Size> channel_idx; ///< indices into the reference-channel vector for this kit's channels
    };

    /// Build a KitLayout for every supported kit (parallel to supportedKits()): its reporter region
    /// (lowest..highest channel m/z, widened by @p buffer on each side) and the ref-channel indices of its channels.
    std::vector<KitLayout> buildKitLayouts(const std::map<uint64_t, Size>& key2ref, double buffer)
    {
      std::vector<KitLayout> layouts;
      for (auto mt : IsobaricKitDetection::supportedKits())
      {
        KitLayout kl;
        double lo = std::numeric_limits<double>::max();
        double hi = std::numeric_limits<double>::lowest();
        for (const auto& c : getKitChannels(mt))
        {
          kl.channel_idx.push_back(key2ref.at(massKey(c.second)));
          lo = std::min(lo, c.second);
          hi = std::max(hi, c.second);
        }
        kl.region_lo = lo - buffer;
        kl.region_hi = hi + buffer;
        layouts.push_back(std::move(kl));
      }
      return layouts;
    }

    /// sum of peak intensities in [lo, hi] of an m/z-sorted spectrum
    double rangeSum(const MSSpectrum& s, double lo, double hi)
    {
      double sum = 0.0;
      for (auto it = s.MZBegin(lo); it != s.end() && it->getMZ() <= hi; ++it) { sum += it->getIntensity(); }
      return sum;
    }

    /// Which MS level carries the reporter ions: 3 if the experiment contains any MS3 spectrum
    /// (SPS-MS3 / MultiNotch TMT workflows quantify reporters in MS3), otherwise 2. When MS3 is present,
    /// the MS2 spectra are ignored entirely.
    Size reporterMsLevel(const PeakMap& exp)
    {
      for (const auto& s : exp) { if (s.getMSLevel() == 3) { return 3; } }
      return 2;
    }

    /// Output of the single measurement pass.
    struct Measurement
    {
      std::vector<IsobaricKitDetection::ChannelStats> channel_stats; ///< parallel to refs
      std::vector<Size> kit_pass;   ///< parallel to kit_layouts: # reporter spectra in which the kit explains >= min_region_coverage of its region intensity
      Size n_signal_spectra = 0;    ///< # reporter (MS2 or MS3) spectra with any signal in the (global) reporter region
    };

    /// Single pass over all reporter spectra (those at @p ms_level): centroid profile spectra on the fly,
    /// match each reference channel, aggregate per-channel statistics, and (per kit) count the spectra whose
    /// reporter-region intensity is dominated by that kit's own channels -- the validity gate.
    /// @p refs must be non-empty & sorted.
    Measurement measure(const PeakMap& exp,
                        const std::vector<IsobaricKitDetection::ChannelRef>& refs,
                        const std::vector<double>& tol,
                        const std::vector<KitLayout>& kit_layouts,
                        Size ms_level,
                        const IsobaricKitDetection::Parameters& params)
    {
      const Size n_ref = refs.size();
      Measurement out;
      out.channel_stats.resize(n_ref);
      out.kit_pass.assign(kit_layouts.size(), 0);
      if (n_ref == 0) { return out; }

      std::vector<std::vector<double>> dppm(n_ref);    // per-channel Δppm samples
      std::vector<std::vector<double>> rel_int(n_ref); // per-channel relative-intensity samples
      std::vector<Size> n_pop(n_ref, 0);               // per-channel population count
      std::vector<double> matched(n_ref, 0.0);         // per-spectrum matched intensity (reused each spectrum)
      const double region_lo = refs.front().mz - 0.2;
      const double region_hi = refs.back().mz + 0.2;

      PeakPickerHiRes picker;
      for (const auto& s : exp)
      {
        if (s.getMSLevel() != ms_level) { continue; }

        const MSSpectrum* ps = &s;
        MSSpectrum tmp;
        if (s.getType(true) == SpectrumSettings::SpectrumType::PROFILE)
        {
          picker.pick(s, tmp);
          ps = &tmp;
        }
        else if (!s.isSorted())
        {
          tmp = s;
          tmp.sortByPosition();
          ps = &tmp;
        }

        // total intensity over the whole reporter region (incl. non-reporter ions)
        const double region_sum = rangeSum(*ps, region_lo, region_hi);
        if (region_sum <= 0.0) { continue; }
        ++out.n_signal_spectra;

        std::fill(matched.begin(), matched.end(), 0.0);
        for (Size i = 0; i < n_ref; ++i)
        {
          Int idx = ps->findNearest(refs[i].mz, tol[i]);
          if (idx < 0) { continue; }
          const double obs = (*ps)[idx].getMZ();
          const double in = (*ps)[idx].getIntensity();
          dppm[i].push_back(Math::getPPM(obs, refs[i].mz));
          rel_int[i].push_back(in / region_sum);
          matched[i] = in;
          ++n_pop[i];
        }

        // validity gate: in this spectrum, do the kit's own channels dominate its own reporter region?
        for (Size k = 0; k < kit_layouts.size(); ++k)
        {
          const KitLayout& kl = kit_layouts[k];
          const double kit_region_sum = rangeSum(*ps, kl.region_lo, kl.region_hi);
          if (kit_region_sum <= 0.0) { continue; }
          double kit_int = 0.0;
          for (Size ri : kl.channel_idx) { kit_int += matched[ri]; }
          if (kit_int >= params.min_region_coverage * kit_region_sum) { ++out.kit_pass[k]; }
        }
      }

      for (Size i = 0; i < n_ref; ++i)
      {
        IsobaricKitDetection::ChannelStats& cs = out.channel_stats[i];
        cs.name = refs[i].name;
        cs.theoretical_mz = refs[i].mz;
        cs.n_found = n_pop[i];
        cs.found_fraction = (out.n_signal_spectra > 0) ? static_cast<double>(n_pop[i]) / out.n_signal_spectra : 0.0;
        if (!dppm[i].empty())
        {
          cs.ppm_offset = medianOf(dppm[i]);
          // sd is undefined (NaN) for a single sample -> report 0 in that case
          cs.ppm_spread = (dppm[i].size() >= 2) ? Math::sd(dppm[i].begin(), dppm[i].end()) : 0.0;
          cs.intensity_share = medianOf(rel_int[i]);
        }
      }
      return out;
    }

    /// Score one candidate kit against the global measurement: assemble its channels (with kit-specific
    /// labels, carrying the already-set global 'ok'/outlier flags), count how many @em detected ('ok')
    /// channels it covers, and compute its (un-normalized) overlap score gated by validity.
    /// @p ref_stats must already be classified. @p detected[i] == (ref channel i is 'ok'). @p key2ref maps
    /// a channel mass-key to its index in @p ref_stats. @p layout supplies the kit's reporter region.
    IsobaricKitDetection::KitResult scoreKit(IsobaricKitDetection::MethodType mt,
                                             const std::vector<IsobaricKitDetection::ChannelStats>& ref_stats,
                                             const std::vector<bool>& detected,
                                             Size detected_count,
                                             const std::map<uint64_t, Size>& key2ref,
                                             const KitLayout& layout,
                                             double region_dominance,
                                             const IsobaricKitDetection::Parameters& params)
    {
      IsobaricKitDetection::KitResult kr;
      kr.type = mt;
      kr.region_dominance = region_dominance;
      kr.is_valid = region_dominance >= params.min_valid_spectra_fraction;
      kr.region_low = layout.region_lo;
      kr.region_high = layout.region_hi;

      const auto chans = getKitChannels(mt);
      kr.num_channels = chans.size();

      std::set<uint64_t> kit_keys;
      for (const auto& c : chans)
      {
        const uint64_t key = massKey(c.second);
        kit_keys.insert(key);
        const Size ri = key2ref.at(key);
        kr.channels.push_back(ref_stats[ri]);     // carries the global classification flags
        kr.channels.back().name = c.first;         // use the kit's own channel label (e.g. "131" vs "131N")
      }

      // diagnostic: fraction of the reporter-region signal captured by the kit's 'ok' (detected) channels
      double clean_signal = 0.0;
      for (const auto& ch : kr.channels) { if (!ch.is_outlier) { clean_signal += ch.intensity_share; } }
      kr.clean_signal_fraction = clean_signal;

      // how many detected ('ok') channels does this kit cover?
      Size covered = 0;
      for (Size i = 0; i < ref_stats.size(); ++i)
      {
        if (detected[i] && kit_keys.count(massKey(ref_stats[i].theoretical_mz))) { ++covered; }
      }
      kr.num_covered = covered;
      kr.num_uncovered = detected_count - covered;
      const double score = IsobaricKitDetection::kitScore(kr.num_channels, detected_count, covered);
      // validity gate: a kit whose channels do not dominate its own reporter region is not a candidate
      kr.score = kr.is_valid ? score : 0.0;
      return kr;
    }

    /// Normalize the kits' raw scores (sum to 1) and sort by descending score, breaking ties towards
    /// the more parsimonious (fewer-channel) kit.
    void normalizeAndSortByScore(std::vector<IsobaricKitDetection::KitResult>& results)
    {
      double score_sum = 0.0;
      for (const auto& r : results) { score_sum += r.score; }
      if (score_sum > 0.0)
      {
        for (auto& r : results) { r.score /= score_sum; }
      }
      std::sort(results.begin(), results.end(),
                [](const IsobaricKitDetection::KitResult& a, const IsobaricKitDetection::KitResult& b)
      {
        if (a.score != b.score) { return a.score > b.score; }
        return a.num_channels < b.num_channels;
      });
    }
  } // anonymous namespace

  std::string IsobaricKitDetection::methodName(MethodType mt)
  {
    return std::string(IsobaricQuantitationMethod::methodDisplayName(mt));
  }

  const std::vector<IsobaricKitDetection::MethodType>& IsobaricKitDetection::supportedKits()
  {
    // Derived directly from the IsobaricQuantitationMethod::MethodType enum so that every concrete
    // isobaric method (TMT *and* iTRAQ) is automatically considered - new methods need no change here.
    // UNKNOWN (the disabled/none sentinel) and the SIZE_OF_METHODTYPE marker are excluded.
    static const std::vector<MethodType> kits = []
    {
      std::vector<MethodType> v;
      for (int i = static_cast<int>(MethodType::UNKNOWN) + 1; i < static_cast<int>(MethodType::SIZE_OF_METHODTYPE); ++i)
      {
        v.push_back(static_cast<MethodType>(i));
      }
      return v;
    }();
    return kits;
  }

  std::vector<IsobaricKitDetection::HierarchyNode> IsobaricKitDetection::buildHierarchy()
  {
    const auto& kits = supportedKits();
    const Size n = kits.size();

    // channel-mass key set per kit
    std::vector<std::set<uint64_t>> keys(n);
    for (Size i = 0; i < n; ++i)
    {
      for (const auto& c : getKitChannels(kits[i])) { keys[i].insert(massKey(c.second)); }
    }

    auto isSubset = [](const std::set<uint64_t>& a, const std::set<uint64_t>& b)
    {
      for (uint64_t k : a) { if (b.find(k) == b.end()) { return false; } }
      return true;
    };

    std::vector<HierarchyNode> nodes(n);
    for (Size i = 0; i < n; ++i) { nodes[i].type = kits[i]; }

    for (Size a = 0; a < n; ++a)
    {
      for (Size b = 0; b < n; ++b)
      {
        if (a == b || keys[a].size() >= keys[b].size()) { continue; }
        if (!isSubset(keys[a], keys[b])) { continue; }
        // b is a (proper) superset of a; keep only direct (minimal) supersets:
        // skip if there is an intermediate kit c with a subset c subset b
        bool minimal = true;
        for (Size c = 0; c < n && minimal; ++c)
        {
          if (c == a || c == b) { continue; }
          if (keys[c].size() > keys[a].size() && keys[c].size() < keys[b].size()
              && isSubset(keys[a], keys[c]) && isSubset(keys[c], keys[b]))
          {
            minimal = false;
          }
        }
        if (minimal)
        {
          nodes[a].parents.push_back(kits[b]);
          nodes[b].children.push_back(kits[a]);
        }
      }
    }
    return nodes;
  }

  std::vector<IsobaricKitDetection::ChannelRef> IsobaricKitDetection::referenceChannels()
  {
    // union of every kit's reporter ions, de-duplicated by exact mass-key (identical channels across
    // kits share a constant -> a single entry). The std::map yields them in ascending m/z order.
    std::map<uint64_t, ChannelRef> ref_map;
    for (MethodType mt : supportedKits())
    {
      for (const auto& c : getKitChannels(mt)) { ref_map.emplace(massKey(c.second), ChannelRef{c.first, c.second}); }
    }
    std::vector<ChannelRef> refs;
    refs.reserve(ref_map.size());
    for (const auto& [mass_key, ref] : ref_map) { refs.push_back(ref); }
    return refs;
  }

  std::vector<double> IsobaricKitDetection::channelTolerances(const std::vector<ChannelRef>& refs, double max_tolerance_ppm)
  {
    // <= max_tolerance_ppm, but never more than half the distance to the channel's nearest neighbour
    // (avoids cross-talk). Sparse neighbourhoods (e.g. iTRAQ, ~1 Th apart) keep the full ppm tolerance;
    // the dense TMTpro N/ND/C/CD quartets (~2.9 mDa) get a tight cap.
    const Size n = refs.size();
    std::vector<double> tol(n);
    for (Size i = 0; i < n; ++i)
    {
      double nn = std::numeric_limits<double>::max(); // distance to nearest neighbour (left/right); stays inf for a lone channel
      if (i > 0)     { nn = std::min(nn, refs[i].mz - refs[i - 1].mz); }
      if (i + 1 < n) { nn = std::min(nn, refs[i + 1].mz - refs[i].mz); }
      tol[i] = std::min(Math::ppmToMass(max_tolerance_ppm, refs[i].mz), 0.5 * nn);
    }
    return tol;
  }

  void IsobaricKitDetection::classifyChannels(std::vector<ChannelStats>& channels, const std::vector<double>& channel_tol_ppm, const Parameters& params)
  {
    // presence/abundance FIRST: the strongest channel sets the reference abundance
    double kit_max_found = 0.0;
    for (const auto& ch : channels) { kit_max_found = std::max(kit_max_found, ch.found_fraction); }

    auto is_underpop = [&](const ChannelStats& ch)
    {
      return kit_max_found > 0.0 && ch.found_fraction < params.underpop_factor * kit_max_found;
    };

    // Robust mass-accuracy baselines from the reliable channels only: well-populated (not under-populated)
    // AND with enough samples (>=2) to trust their Δppm. Weak/absent channels cannot distort the baselines.
    std::vector<double> reliable_sd;   // their ppm_spread   (scatter / precision)
    std::vector<double> reliable_off;  // their ppm_offset   (offset / accuracy)
    for (const auto& ch : channels)
    {
      if (ch.n_found >= 2 && !is_underpop(ch))
      {
        reliable_sd.push_back(ch.ppm_spread);
        reliable_off.push_back(ch.ppm_offset);
      }
    }
    const double med_sd = reliable_sd.empty() ? 0.0 : Math::median(reliable_sd.begin(), reliable_sd.end(), false);
    const double mad_sd = reliable_sd.empty() ? 0.0 : Math::MAD(reliable_sd.begin(), reliable_sd.end(), med_sd);
    const bool use_relative_test = reliable_sd.size() >= params.min_channels_for_mad && mad_sd > 0.0;

    // Consensus instrument-calibration offset shared by all real reporter channels (robust median).
    const double consensus_off = reliable_off.empty() ? 0.0 : Math::median(reliable_off.begin(), reliable_off.end(), false);
    const bool use_offset_test = reliable_off.size() >= 2; // need >= 2 channels for a meaningful consensus

    for (Size k = 0; k < channels.size(); ++k)
    {
      ChannelStats& ch = channels[k];
      ch.is_outlier = false;
      ch.outlier_reason.clear();

      if (ch.n_found == 0)
      { // channel never observed -> not part of the data at all
        ch.is_outlier = true;
        ch.outlier_reason = "missing";
        continue;
      }
      if (is_underpop(ch))
      { // present, but far weaker than the best channel
        ch.is_outlier = true;
        ch.outlier_reason = "underpopulated";
        continue;
      }
      // mass-accuracy ('noisy') tests, judged only on the remaining well-populated channels:
      //  - absolute: scatter approaches the uniform-noise level implied by the +-tol matching window
      //              (scale-free; works even for tiny kits with too few channels for robust statistics)
      //  - relative: scatter is a robust (median/MAD) outlier among the reliable channels
      const double tol_ppm = (k < channel_tol_ppm.size()) ? channel_tol_ppm[k] : 0.0;
      const bool noisy_abs = tol_ppm > 0.0 && ch.ppm_spread > params.noise_sd_frac_of_tol * tol_ppm;
      const bool noisy_rel = use_relative_test && ch.ppm_spread > med_sd + params.ppm_outlier_mad * mad_sd;
      if (noisy_abs || noisy_rel)
      {
        ch.is_outlier = true;
        ch.outlier_reason = "high delta-ppm variance";
        continue;
      }
      // mass-accuracy location: a real reporter channel sits at the shared instrument-calibration offset;
      // a coincidental noise match has a random offset and is flagged here.
      if (use_offset_test && std::abs(ch.ppm_offset - consensus_off) > params.offset_consistency_ppm)
      {
        ch.is_outlier = true;
        ch.outlier_reason = "median deltaPPM inconsistent";
      }
    }
  }

  double IsobaricKitDetection::kitScore(Size num_channels, Size detected_count, Size num_covered)
  {
    // Jaccard overlap of the kit's channel set with the detected-channel set:
    //   1.0 exactly when the kit equals the detected set (best & most parsimonious);
    //   penalised both for missing detected channels (kit too small) and for surplus channels (kit too big).
    const double uni = static_cast<double>(num_channels) + detected_count - num_covered;
    return (uni > 0.0) ? static_cast<double>(num_covered) / uni : 0.0;
  }

  std::vector<IsobaricKitDetection::KitResult> IsobaricKitDetection::detect(const PeakMap& exp, const Parameters& params)
  {
    // 1) reference channels (union of all kits) and 2) their per-channel matching tolerances
    const auto refs = referenceChannels();
    if (refs.empty()) { return {}; }
    const auto tol = channelTolerances(refs, params.max_tolerance_ppm);

    std::map<uint64_t, Size> key2ref; // channel mass-key -> index into refs / ref_stats
    for (Size i = 0; i < refs.size(); ++i) { key2ref[massKey(refs[i].mz)] = i; }
    const auto kit_layouts = buildKitLayouts(key2ref, params.kit_region_buffer); // parallel to supportedKits()

    // reporter ions live in MS3 for SPS-MS3 / MultiNotch workflows; otherwise in MS2. If any MS3 spectrum
    // is present, use only MS3 (and ignore MS2); else use MS2.
    const Size ms_level = reporterMsLevel(exp);

    // 3+4) single pass: quantify the reporter region of every reporter spectrum -> per-channel statistics + per-kit coverage
    Measurement meas = measure(exp, refs, tol, kit_layouts, ms_level, params);
    if (meas.n_signal_spectra == 0)
    {
      OPENMS_LOG_INFO << "IsobaricKitDetection: no MS" << ms_level << " spectra with signal in the reporter region ["
                      << (refs.front().mz - 0.2) << ", " << (refs.back().mz + 0.2)
                      << "] Th were found - cannot detect an isobaric kit." << std::endl;
      return {};
    }
    std::vector<ChannelStats>& ref_stats = meas.channel_stats;

    // 5) classify the whole reference set once -> the 'ok' channels are the DETECTED channels
    std::vector<double> ref_tol_ppm(refs.size());
    for (Size i = 0; i < refs.size(); ++i) { ref_tol_ppm[i] = refs[i].mz > 0.0 ? tol[i] / refs[i].mz * 1e6 : 0.0; }
    classifyChannels(ref_stats, ref_tol_ppm, params);
    std::vector<bool> detected(refs.size());
    for (Size i = 0; i < refs.size(); ++i) { detected[i] = !ref_stats[i].is_outlier; }
    const Size detected_count = static_cast<Size>(std::count(detected.begin(), detected.end(), true));

    // 6) score each candidate kit (with its region-dominance), then normalize to scores and sort
    std::vector<KitResult> results;
    const auto& kits = supportedKits();
    results.reserve(kits.size());
    for (Size k = 0; k < kits.size(); ++k)
    {
      const double region_dominance = static_cast<double>(meas.kit_pass[k]) / meas.n_signal_spectra;
      results.push_back(scoreKit(kits[k], ref_stats, detected, detected_count, key2ref, kit_layouts[k], region_dominance, params));
    }
    normalizeAndSortByScore(results);

    logResults_(results, detected_count, meas.n_signal_spectra, ms_level);
    return results;
  }

  void IsobaricKitDetection::logResults_(const std::vector<KitResult>& results, Size detected_count, Size n_signal_spectra, Size ms_level)
  {
    OPENMS_LOG_INFO << "\n-- Isobaric kit detection --\n"
                    << "Reporter spectra used: MS" << ms_level << "\n"
                    << "MS" << ms_level << " spectra with reporter-region signal: " << n_signal_spectra << "\n"
                    << "Detected ('ok') channels in the data: " << detected_count << "\n" << std::endl;

    for (const auto& kr : results)
    {
      OPENMS_LOG_INFO << methodName(kr.type) << "  (" << kr.num_channels << " channels)"
                      << "   score=" << StringUtils::number(kr.score * 100.0, 1) << "%"
                      << "   owns-region[" << StringUtils::number(kr.region_low, 1) << " - " << StringUtils::number(kr.region_high, 1) << "]="
                      << StringUtils::number(kr.region_dominance * 100.0, 1) << "%" << (kr.is_valid ? "" : " [REJECTED]")
                      << "   clean-signal=" << StringUtils::number(kr.clean_signal_fraction * 100.0, 1) << "%"
                      << "   covers " << kr.num_covered << "/" << detected_count << " detected channels"
                      << (kr.num_uncovered > 0 ? "  [too small]" : "") << "\n";
      OPENMS_LOG_INFO << "    channel    theo.m/z    ppm.offset   ppm.spread   found%   int.share   status\n";
      for (const auto& ch : kr.channels)
      {
        // "missing" (never observed) is shown plainly; other anomalies as "outlier: <reason>"
        std::string status = "ok";
        if (ch.is_outlier) { status = (ch.outlier_reason == "missing") ? ch.outlier_reason : ("outlier: " + ch.outlier_reason); }
        OPENMS_LOG_INFO << "    " << StringUtils::fillRight(std::string(ch.name), ' ', 8)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.theoretical_mz, 5), ' ', 11)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.ppm_offset, 2), ' ', 9)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.ppm_spread, 2), ' ', 9)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.found_fraction * 100.0, 1), ' ', 5)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.intensity_share, 4), ' ', 8)
                        << "  " << status
                        << "\n";
      }
      OPENMS_LOG_INFO << std::endl;
    }

    // most parsimonious *valid* kit that covers all detected channels
    const KitResult* parsimonious = nullptr;
    for (const auto& kr : results)
    {
      if (kr.is_valid && detected_count > 0 && kr.num_uncovered == 0)
      {
        if (parsimonious == nullptr || kr.num_channels < parsimonious->num_channels) { parsimonious = &kr; }
      }
    }

    // results are sorted by score; invalid kits have score 0 and sort last
    const bool any_valid = !results.empty() && results.front().is_valid && results.front().score > 0.0;
    if (any_valid)
    {
      OPENMS_LOG_INFO << "Most likely isobaric kit: " << methodName(results.front().type)
                      << "  (score " << StringUtils::number(results.front().score * 100.0, 1) << "%)" << std::endl;
      if (parsimonious != nullptr)
      {
        OPENMS_LOG_INFO << "Most parsimonious kit covering all detected channels: " << methodName(parsimonious->type)
                        << "  (" << parsimonious->num_channels << " channels)" << std::endl;
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No valid isobaric kit detected (no kit's channels dominate its reporter region) - "
                         "the data does not appear to be isobarically labelled." << std::endl;
    }
  }

} // namespace OpenMS
