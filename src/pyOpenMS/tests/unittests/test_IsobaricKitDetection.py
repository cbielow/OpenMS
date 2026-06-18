import unittest

import pyopenms


class TestIsobaricKitDetection(unittest.TestCase):

    def test_supported_kits(self):
        kits = pyopenms.IsobaricKitDetection.supportedKits()
        # all concrete kits: TMT 6/10/11/16/18/32/35 + iTRAQ 4/8
        self.assertEqual(len(kits), 9)
        MT = pyopenms.IsobaricKitDetection.MethodType
        self.assertIn(MT.TMT_6PLEX, kits)
        self.assertIn(MT.TMT_35PLEX, kits)
        self.assertIn(MT.ITRAQ_4PLEX, kits)
        self.assertNotIn(MT.UNKNOWN, kits)

    def test_method_name(self):
        MT = pyopenms.IsobaricKitDetection.MethodType
        self.assertEqual(pyopenms.IsobaricKitDetection.methodName(MT.TMT_16PLEX), "TMT 16-plex (TMTpro)")
        self.assertEqual(pyopenms.IsobaricKitDetection.methodName(MT.ITRAQ_4PLEX), "iTRAQ 4-plex")

    def test_reference_channels(self):
        refs = pyopenms.IsobaricKitDetection.referenceChannels()
        # 35 TMT channels (TMT35 union) + 8 iTRAQ channels (113..119, 121)
        self.assertEqual(len(refs), 43)
        mzs = [r.mz for r in refs]
        self.assertEqual(mzs, sorted(mzs))  # ascending
        self.assertEqual(len(set(mzs)), len(mzs))  # unique
        names = [r.name for r in refs]
        self.assertIn("113", names)     # iTRAQ 8-plex only
        self.assertIn("135CD", names)   # TMT 35-plex only

    def test_channel_tolerances(self):
        refs = pyopenms.IsobaricKitDetection.referenceChannels()
        tol = pyopenms.IsobaricKitDetection.channelTolerances(refs, 30.0)
        self.assertEqual(len(tol), len(refs))
        by_name = {r.name: (tol[i], r.mz) for i, r in enumerate(refs)}
        # iTRAQ 114: sparse neighbourhood -> keeps the full 30 ppm
        t114, mz114 = by_name["114"]
        self.assertAlmostEqual(t114, 30e-6 * mz114, places=9)
        # TMTpro deuterated 128ND: dense quartet -> capped well below 30 ppm
        t128nd, mz128nd = by_name["128ND"]
        self.assertLess(t128nd, 30e-6 * mz128nd)

    def test_determine_present_channels(self):
        ikd = pyopenms.IsobaricKitDetection
        chans = []
        for pop, n in [(1.0, 10), (0.5, 5), (0.2, 2), (0.0, 0)]:
            cs = ikd.ChannelStats()
            cs.population_fraction = pop
            cs.n_populated = n
            chans.append(cs)
        present = ikd.determinePresentChannels(chans, 0.3)
        self.assertEqual(list(present), [True, True, False, False])

    def test_classify_channels(self):
        ikd = pyopenms.IsobaricKitDetection
        params = ikd.Parameters()

        def make(pop, n, sd):
            cs = ikd.ChannelStats()
            cs.population_fraction = pop
            cs.n_populated = n
            cs.stddev_delta_ppm = sd
            return cs

        chans = [make(1.0, 10, 0.5),   # ok
                 make(1.0, 10, 6.0),   # noisy (absolute: 6 > 0.4 * 12)
                 make(1.0, 10, 0.5),   # ok
                 make(1.0, 10, 0.5),   # ok
                 make(0.0, 0, 0.0),    # missing
                 make(0.05, 1, 0.0)]   # underpopulated
        tol_ppm = [12.0] * len(chans)
        out = ikd.classifyChannels(chans, tol_ppm, params)
        self.assertFalse(out[0].is_outlier)
        self.assertTrue(out[1].is_outlier)
        self.assertEqual(out[1].outlier_reason, "high delta-ppm variance")
        self.assertEqual(out[4].outlier_reason, "missing")
        self.assertEqual(out[5].outlier_reason, "underpopulated")

    def test_kit_score(self):
        ikd = pyopenms.IsobaricKitDetection
        self.assertAlmostEqual(ikd.kitScore(11, 11, 11), 1.0)
        self.assertAlmostEqual(ikd.kitScore(10, 11, 10), 10.0 / 11.0)
        self.assertAlmostEqual(ikd.kitScore(16, 11, 11), 11.0 / 16.0)
        self.assertAlmostEqual(ikd.kitScore(4, 0, 0), 0.0)

    def test_build_hierarchy(self):
        nodes = pyopenms.IsobaricKitDetection.buildHierarchy()
        self.assertEqual(len(nodes), 9)
        MT = pyopenms.IsobaricKitDetection.MethodType
        by_type = {n.type: n for n in nodes}
        # iTRAQ 4-plex is a subset of iTRAQ 8-plex
        self.assertIn(MT.ITRAQ_8PLEX, list(by_type[MT.ITRAQ_4PLEX].parents))
        # TMT 35-plex is the TMT superset -> no parents
        self.assertEqual(len(by_type[MT.TMT_35PLEX].parents), 0)

    def test_detect_tmt11(self):
        MT = pyopenms.IsobaricKitDetection.MethodType
        tmt11 = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                 129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144500]
        exp = pyopenms.MSExperiment()
        for _ in range(10):
            s = pyopenms.MSSpectrum()
            s.setMSLevel(2)
            s.setType(pyopenms.SpectrumSettings.SpectrumType.CENTROID)
            mz = tmt11 + [130.5]            # + non-reporter background in the region
            inten = [1000.0] * len(tmt11) + [500.0]
            s.set_peaks((mz, inten))
            exp.addSpectrum(s)

        results = pyopenms.IsobaricKitDetection.detect(exp)
        self.assertTrue(len(results) > 0)
        self.assertEqual(results[0].type, MT.TMT_11PLEX)
        self.assertEqual(results[0].num_unexplained_present, 0)
        self.assertGreater(results[0].ok_signal_fraction, 0.9)
        # too-small TMT 10-plex ranks below TMT 11-plex
        by_type = {r.type: r for r in results}
        self.assertLess(by_type[MT.TMT_10PLEX].probability, by_type[MT.TMT_11PLEX].probability)
        # per-channel stats are populated
        self.assertEqual(len(results[0].channels), 11)


if __name__ == "__main__":
    unittest.main()
