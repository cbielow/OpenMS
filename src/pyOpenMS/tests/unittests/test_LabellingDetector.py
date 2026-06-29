import unittest

import pyopenms


class TestSILACDetector(unittest.TestCase):

    def test_ms2data_struct(self):
        d = pyopenms.MS2Data()
        d.RT = 1.5
        d.mz = 500.25
        d.charge = 2
        d.index = 7
        self.assertAlmostEqual(d.RT, 1.5)
        self.assertAlmostEqual(d.mz, 500.25)
        self.assertEqual(d.charge, 2)
        self.assertEqual(d.index, 7)

    def test_silac_detector_bindings(self):
        sd = pyopenms.SILACDetector()
        self.assertEqual(len(sd.getZScores()), 4)
        self.assertEqual(len(sd.getPValues()), 4)
        self.assertAlmostEqual(sd.getSignificanceLevel(), 0.05)  # default before running
        self.assertFalse(sd.getIsSILAC())
        # empty input raises an OpenMS exception
        with self.assertRaises(Exception):
            sd.detectSILAC([])


class TestLabellingDetector(unittest.TestCase):

    def _tmt11_experiment(self, n=10):
        tmt11 = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                 129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144500]
        exp = pyopenms.MSExperiment()
        for i in range(n):
            s = pyopenms.MSSpectrum()
            s.setMSLevel(2)
            s.setRT(float(i))
            s.setType(pyopenms.SpectrumSettings.SpectrumType.CENTROID)
            mz = tmt11 + [130.5]                # + a non-reporter background peak in the region
            inten = [1000.0] * len(tmt11) + [500.0]
            s.set_peaks((mz, inten))            # no precursors set -> SILAC step is n/a (gracefully handled)
            exp.addSpectrum(s)
        exp.updateRanges()
        return exp

    def test_detect_tmt11(self):
        MT = pyopenms.IsobaricKitDetection.MethodType
        r = pyopenms.LabellingDetector.detect(self._tmt11_experiment())
        # isobaric detected via the unified class
        self.assertTrue(r.isobaric_detected)
        self.assertEqual(r.isobaric_kit, MT.TMT_11PLEX)
        self.assertTrue(len(r.isobaric_candidates) > 0)
        self.assertFalse(r.isLabelFree())
        # SILAC could not be run (no MS2 precursor scans) -> reported as n/a, not a crash
        self.assertFalse(r.silac_applicable)
        self.assertFalse(r.silac_detected)
        # report() is a readable summary
        rep = pyopenms.LabellingDetector.report(r)
        self.assertIn("Labelling detection", rep)
        self.assertIn("TMT 11-plex", rep)

    def test_label_free(self):
        # a single MS1 spectrum with a stray peak -> no labelling at all
        exp = pyopenms.MSExperiment()
        s = pyopenms.MSSpectrum()
        s.setMSLevel(1)
        s.set_peaks(([500.0], [1000.0]))
        exp.addSpectrum(s)
        r = pyopenms.LabellingDetector.detect(exp)
        self.assertFalse(r.isobaric_detected)
        self.assertTrue(r.isLabelFree())


if __name__ == "__main__":
    unittest.main()
