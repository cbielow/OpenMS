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

    /// Output of the single measurement pass.
    struct Measurement
    {
      std::vector<IsobaricKitDetection::ChannelStats> channel_stats; ///< parallel to refs
      std::vector<Size> kit_pass;   ///< parallel to kit_layouts: # MS2 spectra in which the kit explains >= min_region_coverage of its region intensity
      Size n_ms2_signal = 0;        ///< # MS2 spectra with any signal in the (global) reporter region
    };

    /// Single pass over all MS2 spectra: centroid profile spectra on the fly, match each reference channel,
    /// aggregate per-channel statistics, and (per kit) count the spectra whose reporter-region intensity is
    /// dominated by that kit's own channels -- the per-kit validity gate. @p refs must be non-empty & sorted.
    Measurement measure(const PeakMap& exp,
                        const std::vector<IsobaricKitDetection::ChannelRef>& refs,
                        const std::vector<double>& tol,
                        const std::vector<KitLayout>& kit_layouts,
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
        if (s.getMSLevel() != 2) { continue; }

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
        ++out.n_ms2_signal;

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
        cs.expected_mz = refs[i].mz;
        cs.n_populated = n_pop[i];
        cs.population_fraction = (out.n_ms2_signal > 0) ? static_cast<double>(n_pop[i]) / out.n_ms2_signal : 0.0;
        if (!dppm[i].empty())
        {
          cs.median_delta_ppm = medianOf(dppm[i]);
          // sd is undefined (NaN) for a single sample -> report 0 in that case
          cs.stddev_delta_ppm = (dppm[i].size() >= 2) ? Math::sd(dppm[i].begin(), dppm[i].end()) : 0.0;
          cs.median_rel_intensity = medianOf(rel_int[i]);
        }
      }
      return out;
    }

    /// Score one candidate kit against the global measurement: assemble its channels (with kit-specific
    /// labels), classify them, count how many globally-present channels it explains, and compute its
    /// (un-normalized) overlap score. @p key2ref maps a channel mass-key to its index in @p ref_stats.
    IsobaricKitDetection::KitResult scoreKit(IsobaricKitDetection::MethodType mt,
                                             const std::vector<IsobaricKitDetection::ChannelStats>& ref_stats,
                                             const std::vector<double>& tol_da,
                                             const std::vector<bool>& present,
                                             Size present_count,
                                             const std::map<uint64_t, Size>& key2ref,
                                             double valid_fraction,
                                             const IsobaricKitDetection::Parameters& params)
    {
      IsobaricKitDetection::KitResult kr;
      kr.type = mt;
      kr.valid_fraction = valid_fraction;
      kr.is_valid = valid_fraction >= params.min_valid_spectra_fraction;

      const auto chans = getKitChannels(mt);
      kr.num_channels = chans.size();

      std::set<uint64_t> kit_keys;
      std::vector<double> channel_tol_ppm; // matching tolerance in ppm, parallel to kr.channels
      channel_tol_ppm.reserve(chans.size());
      for (const auto& c : chans)
      {
        const uint64_t key = massKey(c.second);
        kit_keys.insert(key);
        const Size ri = key2ref.at(key);
        kr.channels.push_back(ref_stats[ri]);
        kr.channels.back().name = c.first; // use the kit's own channel label (e.g. "131" vs "131N")
        const double mz = ref_stats[ri].expected_mz;
        channel_tol_ppm.push_back(mz > 0.0 ? tol_da[ri] / mz * 1e6 : 0.0);
      }

      IsobaricKitDetection::classifyChannels(kr.channels, channel_tol_ppm, params);

      // diagnostic: fraction of the reporter-region signal captured by the kit's 'ok' (non-outlier) channels
      double ok_signal = 0.0;
      for (const auto& ch : kr.channels) { if (!ch.is_outlier) { ok_signal += ch.median_rel_intensity; } }
      kr.ok_signal_fraction = ok_signal;

      Size explained = 0;
      for (Size i = 0; i < ref_stats.size(); ++i)
      {
        if (present[i] && kit_keys.count(massKey(ref_stats[i].expected_mz))) { ++explained; }
      }
      kr.num_explained = explained;
      kr.num_unexplained_present = present_count - explained;
      const double score = IsobaricKitDetection::kitScore(kr.num_channels, present_count, explained);
      // validity gate: a kit whose channels do not dominate its own reporter region is not a candidate
      kr.probability = kr.is_valid ? score : 0.0;
      return kr;
    }

    /// Normalize the kits' raw scores into probabilities (sum to 1) and sort by descending probability,
    /// breaking ties towards the more parsimonious (fewer-channel) kit.
    void normalizeAndSortByProbability(std::vector<IsobaricKitDetection::KitResult>& results)
    {
      double score_sum = 0.0;
      for (const auto& r : results) { score_sum += r.probability; }
      if (score_sum > 0.0)
      {
        for (auto& r : results) { r.probability /= score_sum; }
      }
      std::sort(results.begin(), results.end(),
                [](const IsobaricKitDetection::KitResult& a, const IsobaricKitDetection::KitResult& b)
      {
        if (a.probability != b.probability) { return a.probability > b.probability; }
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

  std::vector<bool> IsobaricKitDetection::determinePresentChannels(const std::vector<ChannelStats>& channel_stats, double present_fraction)
  {
    double max_pop = 0.0;
    for (const auto& cs : channel_stats) { max_pop = std::max(max_pop, cs.population_fraction); }

    std::vector<bool> present(channel_stats.size(), false);
    for (Size i = 0; i < channel_stats.size(); ++i)
    {
      const auto& cs = channel_stats[i];
      present[i] = (max_pop > 0.0) && (cs.n_populated > 0) && (cs.population_fraction >= present_fraction * max_pop);
    }
    return present;
  }

  void IsobaricKitDetection::classifyChannels(std::vector<ChannelStats>& channels, const std::vector<double>& channel_tol_ppm, const Parameters& params)
  {
    // presence/abundance FIRST: the strongest channel of the kit sets the reference abundance
    double kit_max_pop = 0.0;
    for (const auto& ch : channels) { kit_max_pop = std::max(kit_max_pop, ch.population_fraction); }

    auto is_underpop = [&](const ChannelStats& ch)
    {
      return kit_max_pop > 0.0 && ch.population_fraction < params.underpop_factor * kit_max_pop;
    };

    // Robust mass-accuracy baselines from the reliable channels only: well-populated (not under-populated)
    // AND with enough samples (>=2) to trust their Δppm. Weak/absent channels cannot distort the baselines.
    std::vector<double> reliable_sd;   // their Δppm std-dev  (scatter / precision)
    std::vector<double> reliable_off;  // their median Δppm   (offset / accuracy)
    for (const auto& ch : channels)
    {
      if (ch.n_populated >= 2 && !is_underpop(ch))
      {
        reliable_sd.push_back(ch.stddev_delta_ppm);
        reliable_off.push_back(ch.median_delta_ppm);
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

      if (ch.n_populated == 0)
      { // channel never observed -> not part of the data at all
        ch.is_outlier = true;
        ch.outlier_reason = "missing";
        continue;
      }
      if (is_underpop(ch))
      { // present, but far weaker than the kit's best channel
        ch.is_outlier = true;
        ch.outlier_reason = "underpopulated";
        continue;
      }
      // mass-accuracy ('noisy') tests, judged only on the remaining well-populated channels:
      //  - absolute: scatter approaches the uniform-noise level implied by the +-tol matching window
      //              (scale-free; works even for tiny kits with too few channels for robust statistics)
      //  - relative: scatter is a robust (median/MAD) outlier among this kit's reliable channels
      const double tol_ppm = (k < channel_tol_ppm.size()) ? channel_tol_ppm[k] : 0.0;
      const bool noisy_abs = tol_ppm > 0.0 && ch.stddev_delta_ppm > params.noise_sd_frac_of_tol * tol_ppm;
      const bool noisy_rel = use_relative_test && ch.stddev_delta_ppm > med_sd + params.ppm_outlier_mad * mad_sd;
      if (noisy_abs || noisy_rel)
      {
        ch.is_outlier = true;
        ch.outlier_reason = "high delta-ppm variance";
        continue;
      }
      // mass-accuracy location: a real reporter channel sits at the shared instrument-calibration offset;
      // a coincidental noise match has a random offset and is flagged here.
      if (use_offset_test && std::abs(ch.median_delta_ppm - consensus_off) > params.offset_consistency_ppm)
      {
        ch.is_outlier = true;
        ch.outlier_reason = "median deltaPPM inconsistent";
      }
    }
  }

  double IsobaricKitDetection::kitScore(Size num_channels, Size present_count, Size num_explained)
  {
    // Jaccard overlap of the kit's channel set with the present-channel set:
    //   1.0 exactly when the kit equals the present set (best & most parsimonious);
    //   penalised both for missing present channels (kit too small) and for surplus channels (kit too big).
    const double uni = static_cast<double>(num_channels) + present_count - num_explained;
    return (uni > 0.0) ? static_cast<double>(num_explained) / uni : 0.0;
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

    // 3+4) single pass: quantify the reporter region of every MS2 spectrum -> per-channel statistics + per-kit coverage
    const Measurement meas = measure(exp, refs, tol, kit_layouts, params);
    if (meas.n_ms2_signal == 0)
    {
      OPENMS_LOG_INFO << "IsobaricKitDetection: no MS2 spectra with signal in the reporter region ["
                      << (refs.front().mz - 0.2) << ", " << (refs.back().mz + 0.2)
                      << "] Th were found - cannot detect an isobaric kit." << std::endl;
      return {};
    }
    const auto& ref_stats = meas.channel_stats;

    // 5) globally 'present' channels
    const auto present = determinePresentChannels(ref_stats, params.present_fraction);
    const Size present_count = static_cast<Size>(std::count(present.begin(), present.end(), true));

    // 6) score each candidate kit (with its validity fraction), then normalize to probabilities and sort
    std::vector<KitResult> results;
    const auto& kits = supportedKits();
    results.reserve(kits.size());
    for (Size k = 0; k < kits.size(); ++k)
    {
      const double valid_fraction = static_cast<double>(meas.kit_pass[k]) / meas.n_ms2_signal;
      results.push_back(scoreKit(kits[k], ref_stats, tol, present, present_count, key2ref, valid_fraction, params));
    }
    normalizeAndSortByProbability(results);

    logResults_(results, present_count, meas.n_ms2_signal);
    return results;
  }

  void IsobaricKitDetection::logResults_(const std::vector<KitResult>& results, Size present_count, Size n_ms2_signal)
  {
    OPENMS_LOG_INFO << "\n-- Isobaric kit detection --\n"
                    << "MS2 spectra with reporter-region signal: " << n_ms2_signal << "\n"
                    << "Channels detected as present in the data: " << present_count << "\n" << std::endl;

    for (const auto& kr : results)
    {
      OPENMS_LOG_INFO << methodName(kr.type) << "  (" << kr.num_channels << " channels)"
                      << "   probability=" << StringUtils::number(kr.probability * 100.0, 1) << "%"
                      << "   valid=" << StringUtils::number(kr.valid_fraction * 100.0, 1) << "%" << (kr.is_valid ? "" : " [INVALID]")
                      << "   ok-signal=" << StringUtils::number(kr.ok_signal_fraction * 100.0, 1) << "%"
                      << "   explains " << kr.num_explained << "/" << present_count << " present channels"
                      << (kr.num_unexplained_present > 0 ? "  [too small]" : "") << "\n";
      OPENMS_LOG_INFO << "    channel    exp.m/z     med.dppm   sd.dppm   pop%    rel.int    flag\n";
      for (const auto& ch : kr.channels)
      {
        // "missing" (never observed) is shown plainly; other anomalies as "outlier: <reason>"
        std::string flag = "ok";
        if (ch.is_outlier) { flag = (ch.outlier_reason == "missing") ? ch.outlier_reason : ("outlier: " + ch.outlier_reason); }
        OPENMS_LOG_INFO << "    " << StringUtils::fillRight(std::string(ch.name), ' ', 8)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.expected_mz, 5), ' ', 11)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.median_delta_ppm, 2), ' ', 9)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.stddev_delta_ppm, 2), ' ', 7)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.population_fraction * 100.0, 1), ' ', 5)
                        << "  " << StringUtils::fillLeft(StringUtils::number(ch.median_rel_intensity, 4), ' ', 8)
                        << "  " << flag
                        << "\n";
      }
      OPENMS_LOG_INFO << std::endl;
    }

    // most parsimonious *valid* kit that explains all present channels
    const KitResult* parsimonious = nullptr;
    for (const auto& kr : results)
    {
      if (kr.is_valid && present_count > 0 && kr.num_unexplained_present == 0)
      {
        if (parsimonious == nullptr || kr.num_channels < parsimonious->num_channels) { parsimonious = &kr; }
      }
    }

    // results are sorted by probability; invalid kits have probability 0 and sort last
    const bool any_valid = !results.empty() && results.front().is_valid && results.front().probability > 0.0;
    if (any_valid)
    {
      OPENMS_LOG_INFO << "Most likely isobaric kit: " << methodName(results.front().type)
                      << "  (probability " << StringUtils::number(results.front().probability * 100.0, 1) << "%)" << std::endl;
      if (parsimonious != nullptr)
      {
        OPENMS_LOG_INFO << "Most parsimonious kit explaining all present channels: " << methodName(parsimonious->type)
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
