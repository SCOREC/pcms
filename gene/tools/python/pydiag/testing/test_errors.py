import unittest

import numpy as np

import pydiag.utils.errors as err


class TestCorrtimeFuncs(unittest.TestCase):

    def test_autocorrtime_nocalc(self):
        """Test the cases where no calculation should occur """
        self.assertEqual(err.autocorrtime([], []), 0)
        self.assertEqual(err.autocorrtime([0.5],[1]), 0)
        self.assertEqual(err.autocorrtime([0.5, 1], [0.5, 1]), 0)
        np.testing.assert_array_equal(err.autocorrtime(np.ones((1, 10)), [1]), np.zeros(10))

    def test_autocorrtime_1d_const(self):
        """ Test the edge case of a constant input """
        times = np.arange(10)
        data = np.ones(10)
        self.assertEqual(err.autocorrtime(data, times), 0)

    def test_autocorrtime_1d_calc(self):
        """ Test an actual computation """
        times = np.arange(1000)
        data = 5*np.arange(1000) + 5
        self.assertEqual(err.autocorrtime(data, times), 217)

    def test_autocorrtime_2d_calc(self):
        times = np.arange(1000)
        data = np.outer(np.arange(1000), np.ones(100))
        ref = 217*np.ones(100)
        np.testing.assert_array_equal(err.autocorrtime(data, times), ref)

    def test_autocorrtime_2d_const(self):
        """ Test the edge case of a constant input """
        times = np.arange(1000)
        data = np.ones((1000, 100))
        ref = np.zeros(100)
        np.testing.assert_array_equal(err.autocorrtime(data, times), ref)

    def test_autocorrtime_3d_const(self):
        """ Test the edge case of a constant input in 3d """
        times = np.arange(1000)
        data = np.ones((1000, 1, 10))
        ref = np.zeros((1, 10))
        np.testing.assert_array_equal(err.autocorrtime(data, times), ref)


class TestWindowErr(unittest.TestCase):

    def test_windowerr_nocalc(self):
        np.testing.assert_equal(err.windowerr([], []), (0, 0.0))
        np.testing.assert_equal(err.windowerr([], [], n_win=10), (0, 0.0))
        np.testing.assert_equal(err.windowerr([0, 1], [0, 1]), (0, 0.0))
        np.testing.assert_equal(err.windowerr(np.ones(10), np.arange(10)), (0, 0.0))

    def test_windowerr_calc(self):
        times = np.arange(1000)
        data = np.ones(1000)
        for i in range(1, 1000):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times), [1.4679E-3, 1])

    def test_windowerr_calc_nd(self):
        times = np.arange(1000)
        data = np.ones((1000, 100))
        for i in range(1, 1000):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times)[0], 1.4679E-3*np.ones(100))

    def test_windowerr_calc_nwin(self):
        times = np.arange(1000)
        data = np.ones(1000)
        for i in range(1, 1000):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times, n_win=10), [1.694E-3, 0])

    def test_windowerr_few(self):
        """ Not enough time steps for 5*tcor windows"""
        times = np.arange(20)
        data = np.ones(20)
        for i in range(1, 20):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times), [8.79297E-2, 1])

    def test_windowerr_toofew(self):
        """ Less than 6 2*tcor windows"""
        times = np.arange(12)
        data = np.ones(12)
        for i in range(1, 12):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times), [0, 1])

    def test_windowerr_calc_toomanywin(self):
        times = np.arange(100)
        data = np.ones(100)
        for i in range(1, 100):
            data[i, ...] = 0.5 * data[i-1, ...]
        data += 2.5
        np.testing.assert_almost_equal(err.windowerr(data, times, n_win=101), [0, 0])

if __name__ == '__main__':
    unittest.main()
