import unittest
import numpy as np

from pydiag.utils.comm import CommonData
from pydiag.data.nrgdata import NrgFile, gluenrgdata


class TestNrgData(unittest.TestCase):
    def setUp(self):
        """ Read a test run"""
        self.cm = CommonData(".dat", -1, -2)
        self.nrg = NrgFile("nrg.dat", self.cm)
        self.nrg.generate_timeseries()

    def test_timefield(self):
        """ Make sure the time stamps are read correctly"""
        times = self.nrg.timearray
        np.testing.assert_almost_equal(times, np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538]))

    def test_nrgcolumn(self):
        """ Test if the file is correctly parsed"""
        np.testing.assert_almost_equal(self.nrg.dataarray[:, 0, 6],
                                       [-2.5089000e-12, 3.5936000e-03, 7.1263000e-03, 1.0543000e-02,
                                        1.3797000e-02, 1.6873000e-02])

    def test_gluenrgdata(self):
        """ Test the coupling routine"""
        cm2 = CommonData(".dat", 0.3, -2)
        nrg2 = NrgFile("nrg.dat", cm2)
        nrg2.generate_timeseries()
        glue = gluenrgdata([self.nrg, nrg2])
        times = glue.timearray
        np.testing.assert_almost_equal(times, np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538]))
        np.testing.assert_almost_equal(glue.dataarray[:, 0, 6],
                                       [-2.5089000e-12, 3.5936000e-03, 7.1263000e-03, 1.0543000e-02,
                                        1.3797000e-02, 1.6873000e-02])
