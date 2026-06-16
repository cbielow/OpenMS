// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  /**
    @brief Single source of truth for the high-accuracy reporter-ion m/z of every TMT/TMTpro channel.

    Every TMT(pro) channel that can occur in any of the supported plex kits (TMT 6/10/11/16/18/32/35-plex)
    is listed here exactly once. The concrete IsobaricQuantitationMethod implementations
    (e.g. TMTSixPlexQuantitationMethod) populate their channel centers from these constants so that the
    same physical reporter ion always has a bit-identical m/z regardless of the kit it appears in.
    This avoids floating-point equality problems when comparing/relating channels across methods (e.g.
    when computing the kit subset/superset hierarchy in TMTPlexDetection).

    @note Values are stored as @c double, not @c float: the masses carry 9 significant digits
          (e.g. 126.127726) which exceeds the ~7 digits a 32-bit @c float can represent. Using @c float
          would silently round the values and break the exact-equality contract above.

    Channel naming follows the Thermo convention: the nominal mass plus an N/C suffix for the
    15N / 13C isobars, and an additional D suffix for the deuterated (TMTpro 32/35-plex) reagents,
    e.g. for nominal mass 128: 128N &lt; 128ND &lt; 128C &lt; 128CD in ascending m/z.

    Reporter-ion masses taken from the Thermo TMT / TMTpro reagent documentation and UniMod.
  */
  namespace TMTMasses
  {
    // classic TMT + TMTpro shared backbone (nominal 126..134)
    static constexpr double TMT_126   = 126.127726;
    static constexpr double TMT_127N  = 127.124761;
    static constexpr double TMT_127C  = 127.131081;
    static constexpr double TMT_128N  = 128.128116;
    static constexpr double TMT_128C  = 128.134436;
    static constexpr double TMT_129N  = 129.131471;
    static constexpr double TMT_129C  = 129.137790;
    static constexpr double TMT_130N  = 130.134825;
    static constexpr double TMT_130C  = 130.141145;
    static constexpr double TMT_131N  = 131.138180;
    static constexpr double TMT_131C  = 131.144500;

    // TMTpro extension (nominal 132..135), used by TMT 16/18/32/35-plex
    static constexpr double TMT_132N  = 132.141535;
    static constexpr double TMT_132C  = 132.147855;
    static constexpr double TMT_133N  = 133.144890;
    static constexpr double TMT_133C  = 133.151210;
    static constexpr double TMT_134N  = 134.148245;
    static constexpr double TMT_134C  = 134.154565; ///< TMT 18-plex
    static constexpr double TMT_135N  = 135.151600; ///< TMT 18-plex

    // deuterated TMTpro reagents (TMT 32/35-plex only)
    static constexpr double TMT_127D  = 127.134003;
    static constexpr double TMT_128ND = 128.131038;
    static constexpr double TMT_128CD = 128.137358;
    static constexpr double TMT_129ND = 129.134393;
    static constexpr double TMT_129CD = 129.140713;
    static constexpr double TMT_130ND = 130.137748;
    static constexpr double TMT_130CD = 130.144068;
    static constexpr double TMT_131ND = 131.141103;
    static constexpr double TMT_131CD = 131.147423;
    static constexpr double TMT_132ND = 132.144458;
    static constexpr double TMT_132CD = 132.150778;
    static constexpr double TMT_133ND = 133.147813;
    static constexpr double TMT_133CD = 133.154133;
    static constexpr double TMT_134ND = 134.151171;
    static constexpr double TMT_134CD = 134.157491;
    static constexpr double TMT_135ND = 135.154526;
    static constexpr double TMT_135CD = 135.160846;
  } // namespace TMTMasses
} // namespace OpenMS
