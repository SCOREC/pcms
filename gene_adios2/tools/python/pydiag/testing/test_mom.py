import unittest
import numpy as np

from pydiag.data.momlib import MomFile
from pydiag.utils.comm import CommonData


class TestMomFile(unittest.TestCase):
    def setUp(self):
        self.filename = "mom_ions.dat"
        self.cm = CommonData(".dat", -1, -2)
        self.testobj = MomFile(self.filename, self.cm)

    def test_timefield(self):
        testtimes = np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538])
        np.testing.assert_almost_equal(self.testobj.timearray, testtimes)

    def test_minmaxtime(self):
        np.testing.assert_almost_equal(self.testobj.get_minmaxtime(), [0, 0.538])

    def test_set_time_wrong(self):
        with np.testing.assert_raises(ValueError):
            self.testobj.set_time(1000)

    def test_set_time(self):
        self.testobj.set_time(self.testobj.timearray[1])
        np.testing.assert_almost_equal(self.testobj.time, 0.1076)

    def test_offset_zero(self):
        self.testobj.set_time(self.testobj.timearray[0])
        np.testing.assert_equal(self.testobj.offset(0), 4 + 8 + 4 + 4)

    def test_offset(self):
        self.testobj.set_time(self.testobj.timearray[1])
        # Go forward by 1 timestep, then take the second array
        # position = tesize + tind*(tesize + leapmom) + var*(entrysize + 2*intsize) + intsize
        # leapmom = nmoms*(entrysize + 2*intsize)
        # entrysize = nx0*nky0*nz0*complexsize
        entrysize = 24*8*16*16
        leapmom = 6*(entrysize + 2*4)
        position = 16 + 1*(16 + leapmom) + 1*(entrysize + 2*4) + 4
        np.testing.assert_equal(self.testobj.offset(1), position)

    def test_readvar(self):
        self.testobj.set_time(self.testobj.timearray[1])
        readmom = self.testobj.readvar(1)
        # The full vector is a bit large to compare, so we resort to the sum and shape
        refsum = -0.0835294011177 - 0.726135371553j
        np.testing.assert_almost_equal(np.sum(readmom), refsum)
        np.testing.assert_equal(readmom.shape, (24, 8, 16))

    def test_tpar(self):
        """ Test the routine to access one of the mom arrays"""
        refsum = -0.0835294011177 - 0.726135371553j
        self.testobj.set_time(self.testobj.timearray[1])
        readmom = self.testobj.tpar()
        # The heat source should produce no momentum
        np.testing.assert_allclose(np.sum(readmom), refsum, atol=1e-15)
        self.testobj.file.close()  # Requesting the same data again should not cause file access
        readmom = self.testobj.tpar()
        np.testing.assert_allclose(np.sum(readmom), refsum, atol=1e-15)
