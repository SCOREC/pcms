import unittest
import numpy as np

from pydiag.data.srcmom_data import SrcmomFile
from pydiag.utils.comm import CommonData


class TestSrcmomFile(unittest.TestCase):

    def setUp(self):
        self.filename = "srcmom_ions.dat"
        self.cm = CommonData(".dat", -1, -2)
        self.testobj = SrcmomFile(self.filename, self.cm)

    def test_timefield(self):
        testtimes = np.array([0.,  0.2152,  0.4304])
        np.testing.assert_almost_equal(self.testobj.timearray, testtimes)

    def test_minmaxtime(self):
        np.testing.assert_almost_equal(self.testobj.get_minmaxtime(), [0, 0.4304])

    def test_set_time_wrong(self):
        with np.testing.assert_raises(ValueError):
            self.testobj.set_time(1000)

    def test_set_time(self):
        self.testobj.set_time(self.testobj.timearray[1])
        np.testing.assert_almost_equal(self.testobj.time, 0.2152)

    def test_offset_zero(self):
        self.testobj.set_time(self.testobj.timearray[0])
        np.testing.assert_equal(self.testobj.offset(0, 0), 4+8+4+4)

    def test_offset(self):
        self.testobj.set_time(self.testobj.timearray[1])
        # Go forward by one timesteps, then take the second array in the second block
        # position = tesize + (intsize + 3*nx0*realsize + intsize)*3 +
        #            tesize + (intsize + 3*nx0*realsize + intsize) + intsize+nx0*realsize
        position = 1*(16 + (4 + 3*24*8 + 4)*3) + 16 + (4 + 3*24*8 + 4) + 4 +24*8
        np.testing.assert_equal(self.testobj.offset(1, 1), position)

    def test_readvar(self):
        self.testobj.set_time(self.testobj.timearray[1])
        readsrc = self.testobj.readvar(1, 0)
        # The heat source should produce no momentum
        np.testing.assert_allclose(readsrc, np.zeros(self.testobj.cm.pnt.nx0), atol=1e-15)

    def test_srcmoment(self):
        self.testobj.set_time(self.testobj.timearray[1])
        readsrc = self.testobj.srcmoment("mom", "ck_heat")
        # The heat source should produce no momentum
        np.testing.assert_allclose(readsrc, np.zeros(self.testobj.cm.pnt.nx0), atol=1e-15)
        #self.testobj.file.close()  # Requesting the same data again should not cause file access
        # Not valid if mmaps are used for the files
        readsrc = self.testobj.srcmoment("mom", "ck_heat")
        np.testing.assert_allclose(readsrc, np.zeros(self.testobj.cm.pnt.nx0), atol=1e-15)
