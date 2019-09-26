import unittest
import numpy as np

from pydiag.data.datafiles import RunDataFiles
from pydiag.data.slices import MomFieldSlice
from pydiag.utils.comm import CommonData, DiagSpace


class TestSlices(unittest.TestCase):

    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.rdf = RunDataFiles(self.cm)

    def test_slice_times(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True)
        sl = MomFieldSlice(self.cm, "tpar", "ions", ds, self.rdf)
        sl.generate_timeseries()
        testtimes = np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538])
        np.testing.assert_almost_equal(sl.timearray, testtimes)

    def test_full3d_slice(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True)
        sl = MomFieldSlice(self.cm, "tpar", "ions", ds, self.rdf)
        sl.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(sl.dataarray), 4.111553396866159)

    def test_full3d_slice_phi(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True)
        sl = MomFieldSlice(self.cm, "phi", "ions", ds, self.rdf)
        sl.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(sl.dataarray), 46.96255522181961)

    def test_3d_slice(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True, xrange=(2, 5),
                       yrange=(1, 9), zrange=(0, -3))
        sl = MomFieldSlice(self.cm, "tpar", "ions", ds, self.rdf)
        sl.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(sl.dataarray), 0.0021708657338445846)

    def test_single3d_slice(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True, xrange=(2, 3),
                       yrange=(1, 2), zrange=(-1, None))
        sl = MomFieldSlice(self.cm, "tpar", "ions", ds, self.rdf)
        sl.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(sl.dataarray), 9.613613409651325e-08)

    def test_1d_slice(self):
        ds = DiagSpace(self.cm.spatialgrid, x_fourier=False, y_fourier=True, xrange=(2, 12),
                       yavg=True, zavg=True)
        sl = MomFieldSlice(self.cm, "tpar", "ions", ds, self.rdf)
        sl.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(sl.dataarray), 0.23217094726463758)
