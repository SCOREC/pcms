import unittest
import numpy as np
import pydiag.utils.fourier as fourier


class TestSpatialFourier(unittest.TestCase):

    def test_kx_to_x(self):
        test = np.ones(10)
        ref = np.zeros(10)
        ref[0] = 10
        np.testing.assert_almost_equal(fourier.kx_to_x(test, 10, axis=0), ref)

    def test_x_to_kx(self):
        test = np.ones(10)
        ref = np.zeros(10)
        ref[0] = 10
        np.testing.assert_almost_equal(fourier.x_to_kx(test, 10, axis=0), ref)

    def test_ky_to_y(self):
        test = np.ones(10)
        ref = np.array(
                [19., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1.,
                 -1., 1.])
        np.testing.assert_almost_equal(fourier.ky_to_y(test, 10, axis=0), ref)

    def test_y_to_ky(self):
        test = np.ones(10)
        ref = np.zeros(6)
        ref[0] = 10
        np.testing.assert_almost_equal(fourier.y_to_ky(test, 10, axis=0), ref)
