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

  std::vector<IsobaricKitDetection::KitResult> IsobaricKitDetection::detect(const PeakMap& exp, const Parameters& params)
  {
    // ---------------------------------------------------------------------
    // 1) reference channels = union of all reporter ions over all kits
    // ---------------------------------------------------------------------
    struct RefChannel { std::string name; double mz; };
    std::map<uint64_t, RefChannel> ref_map; // keyed by mass-key -> unique channel
    for (MethodType mt : supportedKits())
    {
      for (const auto& c : getKitChannels(mt)) { ref_map.emplace(massKey(c.second), RefChannel{c.first, c.second}); }
    }
    std::vector<RefChannel> refs;  ///< all possible channels, sorted by ascending m/z
    for (const auto& [mass_key, ref_channel] : ref_map) { refs.push_back(ref_channel); } // std::map iterates by key == ascending m/z
    const Size n_ref = refs.size(); ///< number of unique channels across all kits
    if (n_ref == 0) { return {}; }

    // ---------------------------------------------------------------------
    // 2) matching tolerance: <= max_tolerance_ppm, but never more than half the
    //    minimum distance between any two reference channels (avoids cross-talk)
    // ---------------------------------------------------------------------
    double min_dist = std::numeric_limits<double>::max();
    for (Size i = 1; i < n_ref; ++i) { min_dist = std::min(min_dist, refs[i].mz - refs[i - 1].mz); }
    const double half_min = 0.5 * min_dist;

    std::vector<double> tol(n_ref);
    for (Size i = 0; i < n_ref; ++i)
    {
      tol[i] = std::min(Math::ppmToMass(params.max_tolerance_ppm, refs[i].mz), half_min);
    }
    const double region_lo = refs.front().mz - 0.2;
    const double region_hi = refs.back().mz + 0.2;

    // ---------------------------------------------------------------------
    // 3) single pass over all MS2 spectra (profile spectra are centroided on the fly)
    // ---------------------------------------------------------------------
    std::vector<std::vector<double>> dppm(n_ref);    // per-channel Δppm samples
    std::vector<std::vector<double>> rel_int(n_ref); // per-channel relative-intensity samples
    std::vector<Size> n_pop(n_ref, 0);               // per-channel population count
    Size n_ms2_signal = 0;                           // MS2 spectra that carried any signal in reporter region (includes non-reporter ions)

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
      double region_sum = 0.0;
      for (auto it = ps->MZBegin(region_lo); it != ps->end() && it->getMZ() <= region_hi; ++it)
      {
        region_sum += it->getIntensity();
      }
      if (region_sum <= 0.0) { continue; }
      ++n_ms2_signal;

      // quantify each candidate channel
      for (Size i = 0; i < n_ref; ++i)
      {
        Int idx = ps->findNearest(refs[i].mz, tol[i]);
        if (idx < 0) { continue; }
        const double obs = (*ps)[idx].getMZ();
        dppm[i].push_back(Math::getPPM(obs, refs[i].mz));
        rel_int[i].push_back((*ps)[idx].getIntensity() / region_sum);
        ++n_pop[i];
      }
    }

    if (n_ms2_signal == 0)
    {
      OPENMS_LOG_INFO << "IsobaricKitDetection: no MS2 spectra with signal in the reporter region ["
                      << region_lo << ", " << region_hi << "] Th were found - cannot detect an isobaric kit." << std::endl;
      return {};
    }

    // ---------------------------------------------------------------------
    // 4) aggregate per-channel statistics
    // ---------------------------------------------------------------------
    std::vector<ChannelStats> ref_stats(n_ref);
    for (Size i = 0; i < n_ref; ++i)
    {
      ChannelStats& cs = ref_stats[i];
      cs.name = refs[i].name;
      cs.expected_mz = refs[i].mz;
      cs.n_populated = n_pop[i];
      cs.population_fraction = static_cast<double>(n_pop[i]) / n_ms2_signal;
      if (!dppm[i].empty())
      {
        cs.median_delta_ppm = medianOf(dppm[i]);
        // sd is undefined (NaN) for a single sample -> report 0 in that case
        cs.stddev_delta_ppm = (dppm[i].size() >= 2) ? Math::sd(dppm[i].begin(), dppm[i].end()) : 0.0;
        cs.median_rel_intensity = medianOf(rel_int[i]);
      }
    }

    // ---------------------------------------------------------------------
    // 5) globally 'present' channels: well populated relative to the best channel
    // ---------------------------------------------------------------------
    double max_pop = 0.0;
    for (Size i = 0; i < n_ref; ++i) { max_pop = std::max(max_pop, ref_stats[i].population_fraction); }

    std::vector<bool> present(n_ref, false);
    Size present_count = 0;
    for (Size i = 0; i < n_ref; ++i)
    {
      if (max_pop > 0.0 && n_pop[i] > 0 && ref_stats[i].population_fraction >= params.present_fraction * max_pop)
      {
        present[i] = true;
        ++present_count;
      }
    }

    std::map<uint64_t, Size> key2ref;
    for (Size i = 0; i < n_ref; ++i) { key2ref[massKey(refs[i].mz)] = i; }

    // ---------------------------------------------------------------------
    // 6) score each candidate kit
    // ---------------------------------------------------------------------
    std::vector<KitResult> results;
    for (MethodType mt : supportedKits())
    {
      KitResult kr;
      kr.type = mt;

      const auto chans = getKitChannels(mt);
      kr.num_channels = chans.size();

      std::vector<Size> kit_ref_idx;
      std::set<uint64_t> kit_keys;
      kit_ref_idx.reserve(chans.size());
      for (const auto& c : chans)
      {
        const uint64_t k = massKey(c.second);
        kit_keys.insert(k);
        const Size ri = key2ref[k];
        kit_ref_idx.push_back(ri);
        kr.channels.push_back(ref_stats[ri]);
        kr.channels.back().name = c.first; // use the kit's own channel label (e.g. "131" vs "131N")
      }

      // outlier detection within the kit
      double kit_max_pop = 0.0;
      for (Size ri : kit_ref_idx) { kit_max_pop = std::max(kit_max_pop, ref_stats[ri].population_fraction); }

      std::vector<double> stddevs; // Δppm std-dev of the kit's populated channels
      for (Size k = 0; k < kit_ref_idx.size(); ++k)
      {
        if (present[kit_ref_idx[k]]) { stddevs.push_back(kr.channels[k].stddev_delta_ppm); }
      }
      const double med_sd = stddevs.empty() ? 0.0 : Math::median(stddevs.begin(), stddevs.end(), false);
      const double mad_sd = stddevs.empty() ? 0.0 : Math::MAD(stddevs.begin(), stddevs.end(), med_sd);

      for (Size k = 0; k < kr.channels.size(); ++k)
      {
        ChannelStats& ch = kr.channels[k];
        const Size ri = kit_ref_idx[k];
        if (ch.n_populated == 0)
        { // channel never observed -> not part of the data at all
          ch.is_outlier = true;
          ch.outlier_reason = "missing";
          continue;
        }
        const bool underpop = kit_max_pop > 0.0 && ch.population_fraction < params.underpop_factor * kit_max_pop;
        const bool noisy = present[ri] && mad_sd > 0.0 && ch.stddev_delta_ppm > med_sd + params.ppm_outlier_mad * mad_sd;
        ch.is_outlier = underpop || noisy;
        if (underpop) { ch.outlier_reason = "underpopulated"; }
        else if (noisy) { ch.outlier_reason = "high delta-ppm variance"; }
      }

      // how many globally present channels does this kit explain?
      Size explained = 0;
      for (Size i = 0; i < n_ref; ++i)
      {
        if (present[i] && kit_keys.count(massKey(refs[i].mz))) { ++explained; }
      }
      kr.num_explained = explained;
      kr.num_unexplained_present = present_count - explained;

      // Jaccard similarity between the kit's channel set and the present-channel set:
      //   1.0 exactly when the kit equals the present set (best & most parsimonious)
      //   penalised both for missing present channels (kit too small) and for extra channels (kit too big)
      const double uni = static_cast<double>(kr.num_channels) + present_count - explained;
      kr.probability = (uni > 0.0) ? static_cast<double>(explained) / uni : 0.0;

      results.push_back(std::move(kr));
    }

    // normalize scores into probabilities (sum to 1)
    double score_sum = 0.0;
    for (const auto& r : results) { score_sum += r.probability; }
    if (score_sum > 0.0)
    {
      for (auto& r : results) { r.probability /= score_sum; }
    }

    // sort by probability (desc); break ties towards the more parsimonious (fewer-channel) kit
    std::sort(results.begin(), results.end(), [](const KitResult& a, const KitResult& b)
    {
      if (a.probability != b.probability) { return a.probability > b.probability; }
      return a.num_channels < b.num_channels;
    });

    logResults_(results, present_count, n_ms2_signal);
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

    // most parsimonious kit = smallest kit that explains all present channels
    const KitResult* parsimonious = nullptr;
    for (const auto& kr : results)
    {
      if (present_count > 0 && kr.num_unexplained_present == 0)
      {
        if (parsimonious == nullptr || kr.num_channels < parsimonious->num_channels) { parsimonious = &kr; }
      }
    }

    if (!results.empty())
    {
      OPENMS_LOG_INFO << "Most likely isobaric kit: " << methodName(results.front().type)
                      << "  (probability " << StringUtils::number(results.front().probability * 100.0, 1) << "%)" << std::endl;
    }
    if (parsimonious != nullptr)
    {
      OPENMS_LOG_INFO << "Most parsimonious kit explaining all present channels: " << methodName(parsimonious->type)
                      << "  (" << parsimonious->num_channels << " channels)" << std::endl;
    }
    else if (present_count == 0)
    {
      OPENMS_LOG_INFO << "No isobaric reporter channels were reliably detected." << std::endl;
    }
  }

} // namespace OpenMS
