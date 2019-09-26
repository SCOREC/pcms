import unittest
import numpy as np
import os

from pydiag.utils.comm import CommonData
from pydiag.data.fsamom_data import FSAmomData, FSAmomFile


# TODO: Test for initial profiles, nustar
class TestFSAData(unittest.TestCase):
    def setUp(self):
        """ Read a test run"""
        self.cm = CommonData(".dat", -1, -2)
        self.fsa = FSAmomData(self.cm)

    def test_timefield(self):
        """ Make sure the time stamps are read correctly"""
        self.fsa.getfsadata()
        times = self.fsa.fsaspec[0].timearray
        np.testing.assert_almost_equal(times, np.array([0., 0.2152, 0.4304]))


class TestFSAmomFile(unittest.TestCase):
    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.fsf = FSAmomFile("fsamome_ions.dat", self.cm)

    def test_invalidfile(self):
        with np.testing.assert_raises(IOError):
            temp = FSAmomFile("foobar", self.cm)

    def test_timefield(self):
        times = self.fsf.timearray
        np.testing.assert_almost_equal(times, [0., 0.2152, 0.4304])

    def test_getminmaxtime(self):
        np.testing.assert_almost_equal(self.fsf.get_minmaxtime(), [0, 0.4304])

    def test_singletime(self):
        tempcm = CommonData(".dat", 0.3, 0.3)
        tempfsf = FSAmomFile("fsamome_ions.dat", tempcm)
        tempfsf.generate_timeseries()
        # The class should round to the nearest existing timestamp
        np.testing.assert_almost_equal(tempfsf.timearray, [0.4304])

    def test_changetolasttime(self):
        """ Set the class to read only the last time point, then reset"""
        self.fsf.starttime = -2
        self.fsf.endtime = -2
        self.fsf.generate_timeseries()
        np.testing.assert_almost_equal(self.fsf.timearray, [0.4304])
        self.fsf.starttime = -1
        self.fsf.endtime = -2
        self.fsf.generate_timeseries()

    def test_fsa_profile(self):
        temp = np.array(self.fsf.dataarray)[..., 2]
        np.testing.assert_almost_equal(temp, np.array([
            [-4.55775200e-10, -1.16707300e-09, -2.00889400e-09, -3.00564100e-09, -4.12161300e-09,
             -5.20939000e-09, -6.01479000e-09, -6.25958400e-09, -5.79481800e-09, -4.72287300e-09,
             -3.36452600e-09, -2.08139100e-09, -1.10676600e-09, -4.97818600e-10, -1.86050200e-10,
             -5.69762700e-11, -1.36356400e-11, -1.05629000e-12, 2.71281500e-12, 4.17834700e-12,
             4.97034900e-12, 5.39635300e-12, 6.95603800e-12, 1.61955700e-11],
            [4.30433600e-07, 3.67001800e-07, 2.26327100e-07, -3.04187600e-06, -1.48487400e-06,
                1.43386200e-04, 1.03744600e-03, 3.78266900e-03, 8.25737800e-03, 1.06384300e-02,
                6.08910600e-03, -3.48013500e-03, -9.73957100e-03, -8.67513300e-03, -4.28534900e-03,
                -1.20386000e-03, -1.48271900e-04, 1.17795200e-05, 4.79207200e-06, -6.18159500e-07,
                -4.96055000e-07, -1.91160300e-07, -2.17427700e-08, 2.38119400e-06],
            [3.24394200e-07, 2.22078100e-07, 1.33760600e-07, -3.69611100e-06, 5.10401200e-06,
                2.83310400e-04, 1.98412800e-03, 7.30799200e-03, 1.61989800e-02, 2.11581500e-02,
                1.22904000e-02, -6.80064800e-03, -1.93290700e-02, -1.71939700e-02, -8.46640300e-03,
                -2.38131700e-03, -3.06204000e-04, 1.33851300e-05, 5.69722100e-06, -2.05928900e-06,
                -1.06428000e-06, -3.64131900e-07, -1.16575000e-07, 4.79682000e-06]]))
