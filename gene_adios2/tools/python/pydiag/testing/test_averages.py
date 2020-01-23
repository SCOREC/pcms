import unittest

import numpy as np

import pydiag.utils.averages as av
from pydiag.utils.geom import Geometry
from pydiag.utils.comm import CommonData
from pydiag.data.fieldlib import FieldFile


class TestTrapzFuncs(unittest.TestCase):

    def test_mytrapz_singlepoint(self):
        self.assertEqual(av.mytrapz([0.5], [1]), 0.5)

    def test_mytrapz_multipoint(self):
        times = [1, 2, 3]
        data = [1, 2, 3]
        self.assertAlmostEqual(av.mytrapz(data, times), 2.0)

    def test_mytrapz_zeropoint(self):
        self.assertRaises(IndexError, av.mytrapz, [], [])

    def test_mytrapz_length(self):
        test = np.ones((5, 4))
        times = np.ones(4)
        self.assertRaises(ValueError, av.mytrapz, test, times)


class TestSpatialAverage(unittest.TestCase):

    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.geom = Geometry(self.cm)

    def test_zav_ones(self):
        testvar = np.ones((24, 8, 16))
        res = av.z_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(res, np.ones((24, 8)))

    def test_zav_series(self):
        testvar = np.linspace(0, 128, num=3072).reshape((24, 8, 16))
        res = av.z_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(np.sum(res), 12289.224570376951)

    def test_zav_generun(self):
        ff = FieldFile("field.dat", self.cm)
        ff.set_approximate_time(0.538)
        res = av.x_av3d(np.abs(ff.phi())**2, self.geom)[1,6]
        np.testing.assert_almost_equal(res, 0.029694647412238581)

    def test_xav_ones_xyfourier(self):
        geom = Geometry(self.cm)
        geom.cm.x_local = True
        geom.cm.y_local = True
        testvar = np.ones((24, 8, 16))
        res = av.x_av3d(testvar, geom)
        np.testing.assert_almost_equal(res, 24*np.ones((8, 16)))

    def test_xav_ones_xfourier(self):
        geom = Geometry(self.cm)
        geom.cm.x_local = True
        geom.cm.y_local = False
        testvar = np.ones((24, 8, 16))
        res = av.x_av3d(testvar, geom)
        np.testing.assert_almost_equal(res, 47*np.ones((8, 16)))

    def test_xav_ones_nofourier(self):
        geom = Geometry(self.cm)
        geom.cm.x_local = False
        geom.cm.y_local = False
        testvar = np.ones((24, 8, 16))
        res = av.x_av3d(testvar, geom)
        np.testing.assert_almost_equal(res, np.ones((8, 16)))

    def test_yav_ones_nofourier(self):
        geom = Geometry(self.cm)
        geom.cm.x_local = False
        geom.cm.y_local = False
        testvar = np.ones((24, 8, 16))
        res = av.y_av3d(testvar, geom)
        np.testing.assert_almost_equal(res, np.ones((24, 16)))

    def test_yav_ones_fourier(self):
        geom = Geometry(self.cm)
        geom.cm.x_local = False
        geom.cm.y_local = True
        testvar = np.ones((24, 8, 16))
        res = av.y_av3d(testvar, geom)
        np.testing.assert_almost_equal(res, 15*np.ones((24, 16)))

    def test_xzav_ones(self):
        testvar = np.ones((24, 8, 16))
        res = av.xz_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(res, np.ones(8))

    def test_yzav_ones(self):
        testvar = np.ones((24, 8, 16))
        res = av.yz_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(res, 15*np.ones(24))

    def test_xyav_ones(self):
        testvar = np.ones((24, 8, 16))
        res = av.xy_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(res, 15*np.ones(16))

    def test_xyzav_ones(self):
        testvar = np.ones((24, 8, 16))
        res = av.xyz_av3d(testvar, self.geom)
        np.testing.assert_almost_equal(res, 15)