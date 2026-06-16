// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Detects which isobaric labelling kit (TMT/TMTpro or iTRAQ plex) was used in an LC-MS/MS experiment.

    Given an MSExperiment with MS2 spectra, the algorithm quantifies the reporter-ion region
    of each MS2 spectrum against the union of all reporter ions of all supported isobaric kits
    (every concrete IsobaricQuantitationMethod::MethodType, i.e. TMT 6/10/11/16/18/32/35-plex and
    iTRAQ 4/8-plex) and decides, per kit, how well the observed channel pattern matches that kit. It
    returns a probability for each kit (sorted, highest first) and identifies the most parsimonious
    kit (the one with the fewest channels that still explains all channels actually present in the data).

    The supported kits form a subset/superset @em hierarchy (e.g. TMT 16-plex is a subset of
    TMT 32-plex, iTRAQ 4-plex a subset of iTRAQ 8-plex), which is computed at run time from the shared,
    high-accuracy channel masses (see buildHierarchy()). The hierarchy makes the parsimony decision
    robust: the chosen kit is the smallest kit whose channel set covers every present channel.

    Per-channel statistics (median Δppm, Δppm std-dev, population fraction, relative intensity and
    outlier flags) are logged via OPENMS_LOG_INFO for every candidate kit.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI IsobaricKitDetection
  {
  public:
    using MethodType = IsobaricQuantitationMethod::MethodType;

    /// Tunable thresholds for detection. The defaults work for typical Orbitrap reporter-ion data.
    struct OPENMS_DLLAPI Parameters
    {
      /// upper bound on the matching tolerance in ppm (also capped at half the minimum channel distance of the kit)
      double max_tolerance_ppm = 30.0;
      /// a channel counts as 'present' if its population fraction is >= present_fraction * (max population fraction of best channel)
      double present_fraction = 0.3;
      /// within a kit, a channel is an 'underpopulated' outlier if its population fraction < underpop_factor * (max population fraction within that kit)
      double underpop_factor = 0.3;
      /// within a kit, a channel is a 'noisy' outlier if its Δppm std-dev > median(std-dev) + ppm_outlier_mad * MAD(std-dev) over the kit's populated channels
      double ppm_outlier_mad = 3.0;
    };

    /// Per-channel detection statistics, aggregated across all MS2 spectra.
    struct OPENMS_DLLAPI ChannelStats
    {
      std::string name;                   ///< channel label, e.g. "129C" (kit-specific, not derivable from an enum)
      double expected_mz = 0.0;           ///< theoretical reporter-ion m/z
      double median_delta_ppm = 0.0;      ///< median of (observed - expected) in ppm, over populated spectra
      double stddev_delta_ppm = 0.0;      ///< std-dev of the Δppm values, over populated spectra
      double median_rel_intensity = 0.0;  ///< median of (channel intensity / total reporter-region intensity), over populated spectra
      double population_fraction = 0.0;   ///< fraction of MS2 spectra in which the channel was found, in [0, 1]
      Size   n_populated = 0;             ///< number of MS2 spectra in which the channel was found
      bool   is_outlier = false;          ///< flagged as an outlier within its kit
      std::string outlier_reason;         ///< human-readable reason if @p is_outlier, else empty
    };

    /// Result for one candidate isobaric kit.
    struct OPENMS_DLLAPI KitResult
    {
      MethodType type = MethodType::UNKNOWN; ///< the kit's MethodType (use methodName() for a display name)
      double probability = 0.0;              ///< score normalized across all kits, in [0, 1] (higher = more likely)
      Size num_channels = 0;                 ///< number of channels this kit defines
      Size num_explained = 0;                ///< number of present channels that this kit contains
      Size num_unexplained_present = 0;      ///< number of present channels NOT in this kit (i.e. kit is too small)
      std::vector<ChannelStats> channels;    ///< per-channel statistics for this kit's channels
    };

    /// A node of the kit subset/superset hierarchy.
    struct OPENMS_DLLAPI HierarchyNode
    {
      MethodType type;                       ///< this kit
      std::vector<MethodType> parents;       ///< direct (minimal) supersets of this kit
      std::vector<MethodType> children;      ///< direct (maximal) subsets of this kit
    };

    /// Human-readable display name of a kit (thin wrapper around IsobaricQuantitationMethod::methodDisplayName()).
    static std::string methodName(MethodType mt);

    /**
      @brief Detect the isobaric kit used in @p exp.

      Iterates over all MS2 spectra (profile spectra are centroided on the fly; @p exp itself is not
      modified), quantifies the reporter-ion region and returns one KitResult per supported kit,
      sorted by @c probability (descending). The most likely / most parsimonious kit is the first
      element. Per-channel statistics are additionally logged via OPENMS_LOG_INFO.

      @param exp Input experiment; only MS2 spectra are used.
      @param params Detection thresholds.
      @return One KitResult per supported kit, sorted by descending probability. Empty if @p exp has no MS2 spectra.
    */
    static std::vector<KitResult> detect(const PeakMap& exp, const Parameters& params = Parameters());

    /// Build the subset/superset hierarchy of all supported isobaric kits from their channel masses.
    /// Currently, this spells out as: TMT6 ⊂ TMT10 ⊂ TMT11 ⊂ TMT16 ⊂ {TMT18, TMT32} ⊂ TMT35, and iTRAQ4 ⊂ iTRAQ8.
    static std::vector<HierarchyNode> buildHierarchy();

    /// The isobaric kits considered by detect()/buildHierarchy(): all concrete IsobaricQuantitationMethod::MethodType values.
    static const std::vector<MethodType>& supportedKits();

  private:
    /// Log per-kit / per-channel statistics and the final decision via OPENMS_LOG_INFO.
    static void logResults_(const std::vector<KitResult>& results, Size present_count, Size n_ms2_signal);
  };

} // namespace OpenMS
