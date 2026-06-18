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

    @par How a channel is classified (per candidate kit)
    Presence/abundance is judged first; mass accuracy only on the channels that survive. A channel is
    "ok" only if it reaches the bottom of this tree:
    @code
      n_populated == 0 ?
       |- yes --------------------------------------------------------> OUTLIER: "missing"
       '- no
           under-populated?  (population_fraction < underpop_factor * kit_max_pop)
            |- yes ---------------------------------------------------> OUTLIER: "underpopulated"
            '- no   (channel is well-populated)
                noisy_abs = stddev_dppm > noise_sd_frac_of_tol * tol_ppm           (absolute: vs the
                                                                                    uniform-noise level
                                                                                    of the +-tol window)
                noisy_rel = use_relative_test                                      (relative: robust
                            && stddev_dppm > median_sd + ppm_outlier_mad * MAD_sd   median/MAD outlier
                                                                                    among reliable chans)
                noisy_abs OR noisy_rel ?
                 |- yes --------------------------------------------> OUTLIER: "high delta-ppm variance"
                 '- no
                     |median_dppm - consensus_offset| > offset_consistency_ppm ?  (all real reporter
                                                                                    channels share one
                                                                                    instrument calibration
                                                                                    offset)
                      |- yes -------------------------------------> OUTLIER: "median deltaPPM inconsistent"
                      '- no --------------------------------------> OK
    @endcode
    where
      - @c kit_max_pop = max population_fraction over the kit's channels;
      - the @em reliable channels (which define @c median_sd, @c MAD_sd and @c consensus_offset) are those
        with @c n_populated @c >= @c 2 AND not under-populated -- weak/absent channels never distort the baseline;
      - @c use_relative_test = (#reliable @c >= @c min_channels_for_mad) AND @c MAD_sd @c > @c 0
        (below that the robust baseline is too unstable, so only the absolute test is used);
      - @c consensus_offset = median of the reliable channels' @c median_delta_ppm (the test needs @c >= @c 2
        reliable channels);
      - @c tol_ppm = per-channel match tolerance in ppm = min(max_tolerance_ppm, half the distance to
        the channel's nearest neighbour among all reference channels). Sparse neighbourhoods (e.g. iTRAQ)
        keep the full ppm tolerance; the dense TMTpro N/ND/C/CD quartets get a tight cap.

    @par Kit validity gate (per-kit, not per-channel)
    Before a kit is ranked, it must be @em valid: in at least @c min_valid_spectra_fraction of the MS2
    spectra, the intensity matched to the kit's own channels must be at least @c min_region_coverage of the
    total intensity in the kit's reporter region (the kit's lowest-to-highest channel m/z, widened by
    @c kit_region_buffer on each side). Kits that fail this (e.g. label-free data whose 126-131 region is
    just peptide-fragment noise) get probability 0 and are never reported as the detected kit. See KitResult::is_valid.

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
      /// absolute mass-accuracy test: a channel is 'noisy' if its Δppm std-dev exceeds noise_sd_frac_of_tol * (matching tolerance in ppm).
      /// Random matches within the +-tol window are ~uniform, with std-dev ~tol/sqrt(3) ~= 0.58*tol, so a real channel sits well below this.
      double noise_sd_frac_of_tol = 0.4;
      /// relative mass-accuracy test: a (reliable) channel is also 'noisy' if its Δppm std-dev > median(std-dev) + ppm_outlier_mad * MAD(std-dev),
      /// where median/MAD are taken over the kit's reliable channels only (see min_channels_for_mad)
      double ppm_outlier_mad = 3.0;
      /// minimum number of reliable (well-populated) channels required before the relative (median/MAD) test is applied;
      /// below this the robust baseline is too unstable and only the absolute test is used
      Size min_channels_for_mad = 4;
      /// median-offset consistency: a channel is an outlier if its median Δppm deviates from the reliable channels'
      /// consensus (their median of median Δppm) by more than this many ppm. Real reporter channels share one
      /// instrument calibration offset, so coincidental noise matches (with random offsets) are caught here.
      double offset_consistency_ppm = 5.0;
      /// kit validity gate (per spectrum): the kit's channels must capture at least this fraction of the total
      /// intensity in the kit's own reporter region for the spectrum to count as 'explained'
      double min_region_coverage = 0.5;
      /// kit validity gate (over spectra): a kit is only considered valid if at least this fraction of MS2 spectra
      /// are 'explained' (see min_region_coverage). Invalid kits get probability 0
      double min_valid_spectra_fraction = 0.5;
      /// outward buffer (in Th) added on each side of a kit's lowest..highest channel m/z to define its reporter region
      double kit_region_buffer = 0.1;
    };

    /// A reporter-ion reference channel: its label and theoretical m/z.
    struct OPENMS_DLLAPI ChannelRef
    {
      std::string name; ///< channel label, e.g. "129C"
      double mz = 0.0;  ///< theoretical reporter-ion m/z
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
      double ok_signal_fraction = 0.0;       ///< fraction [0,1] of the reporter-region signal captured by this kit's non-outlier ('ok') channels (sum of their median relative intensities); a diagnostic only -- NOT used for scoring
      double valid_fraction = 0.0;           ///< fraction [0,1] of MS2 spectra in which this kit's channels capture >= min_region_coverage of the kit's reporter-region intensity (see the validity gate)
      bool is_valid = false;                 ///< whether @p valid_fraction >= min_valid_spectra_fraction; invalid kits get probability 0 and are never reported as detected
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

    /// @name Composable building blocks of detect()
    /// These are the individually-testable, side-effect-free steps that detect() orchestrates.
    /// @{

    /// The union of all reporter ions of all supported kits, de-duplicated and sorted by ascending m/z.
    static std::vector<ChannelRef> referenceChannels();

    /// Per-channel matching tolerance in Th (parallel to @p refs): min(@p max_tolerance_ppm, half the
    /// distance to the channel's nearest neighbour). @p refs must be sorted by ascending m/z.
    static std::vector<double> channelTolerances(const std::vector<ChannelRef>& refs, double max_tolerance_ppm);

    /// Boolean mask (parallel to @p channel_stats) of channels considered 'present' in the data:
    /// populated and with population fraction >= @p present_fraction * (max population fraction over all channels).
    static std::vector<bool> determinePresentChannels(const std::vector<ChannelStats>& channel_stats, double present_fraction);

    /**
      @brief Classify each channel of a single kit as ok / missing / underpopulated / noisy / offset-inconsistent.

      Sets @c is_outlier and @c outlier_reason on every entry of @p channels following the decision tree
      documented on the class. Presence/abundance is decided first; the mass-accuracy tests (Δppm scatter,
      then median-Δppm consistency against the reliable channels' consensus) are applied only to the
      surviving well-populated channels.

      @param[in,out] channels per-channel stats of ONE kit; @c n_populated, @c population_fraction and
                              @c stddev_delta_ppm must be filled in. Only the outlier fields are written.
      @param channel_tol_ppm matching tolerance in ppm for each channel (parallel to @p channels).
      @param params classification thresholds.
    */
    static void classifyChannels(std::vector<ChannelStats>& channels, const std::vector<double>& channel_tol_ppm, const Parameters& params);

    /// Overlap score of a kit's channel set against the present-channel set (Jaccard): num_explained / (num_channels + present_count - num_explained).
    /// 1.0 exactly when the kit equals the present set; penalised for both missing and surplus channels.
    static double kitScore(Size num_channels, Size present_count, Size num_explained);

    /// @}

  private:
    /// Log per-kit / per-channel statistics and the final decision via OPENMS_LOG_INFO.
    static void logResults_(const std::vector<KitResult>& results, Size present_count, Size n_ms2_signal);
  };

} // namespace OpenMS
