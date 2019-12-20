import unittest
import numpy as np

from pydiag.utils.comm import CommonData, SpatialGrid, DiagSpace


class TestCommon(unittest.TestCase):
    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)

    def test_local(self):
        self.assertFalse(self.cm.x_local)
        self.assertTrue(self.cm.y_local)

    def test_nonlinear(self):
        self.assertTrue(self.cm.nonlinear)


class TestSpatialGrid(unittest.TestCase):
    def setUp(self):
        cm = CommonData(".dat", -1, -2)
        self.grid_xgyl = cm.spatialgrid
        # Create mock grids from the parameters to test the other global/local combinations
        self.grid_xgyg = SpatialGrid(cm.pnt, x_local=False, y_local=False)
        self.grid_xlyg = SpatialGrid(cm.pnt, x_local=True, y_local=False)
        self.grid_xlyl = SpatialGrid(cm.pnt, x_local=True, y_local=True)

    def test_xgrid_xglobal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xgyl.x_a, np.array(
                [0.052, 0.09095652, 0.12991304, 0.16886957, 0.20782609, 0.24678261, 0.28573913,
                 0.32469565, 0.36365217, 0.4026087, 0.44156522, 0.48052174, 0.51947826, 0.55843478,
                 0.5973913, 0.63634783, 0.67530435, 0.71426087, 0.75321739, 0.79217391, 0.83113043,
                 0.87008696, 0.90904348, 0.948]))

    def test_xgrid_xglobal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xgyg.x_a, np.array(
                [0.052, 0.09095652, 0.12991304, 0.16886957, 0.20782609, 0.24678261, 0.28573913,
                 0.32469565, 0.36365217, 0.4026087, 0.44156522, 0.48052174, 0.51947826, 0.55843478,
                 0.5973913, 0.63634783, 0.67530435, 0.71426087, 0.75321739, 0.79217391, 0.83113043,
                 0.87008696, 0.90904348, 0.948]))

    def test_xgrid_xlocal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xlyg.x, np.array(
                [-80., -73.33333333, -66.66666667, -60., -53.33333333, -46.66666667, -40.,
                 -33.33333333, -26.66666667, -20., -13.33333333, -6.66666667, 0., 6.66666667,
                 13.33333333, 20., 26.66666667, 33.33333333, 40., 46.66666667, 53.33333333, 60.,
                 66.66666667, 73.33333333]))

    def test_xgrid_xlocal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xlyl.x, np.array(
                [-80., -73.33333333, -66.66666667, -60., -53.33333333, -46.66666667, -40.,
                 -33.33333333, -26.66666667, -20., -13.33333333, -6.66666667, 0., 6.66666667,
                 13.33333333, 20., 26.66666667, 33.33333333, 40., 46.66666667, 53.33333333, 60.,
                 66.66666667, 73.33333333]))

    def test_ygrid_xglobal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xgyl.y, np.array(
                [-6.6742992, -5.8400118, -5.0057244, -4.171437, -3.3371496, -2.5028622, -1.6685748,
                 -0.8342874, 0., 0.8342874, 1.6685748, 2.5028622, 3.3371496, 4.171437, 5.0057244,
                 5.8400118]))

    def test_ygrid_xglobal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xgyg.y, np.array(
                [-6.6742992, -5.0057244, -3.3371496, -1.6685748, 0., 1.6685748, 3.3371496,
                 5.0057244]))

    def test_ygrid_xlocal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xlyg.y, np.array(
                [-6.6742992, -5.0057244, -3.3371496, -1.6685748, 0., 1.6685748, 3.3371496,
                 5.0057244]))

    def test_ygrid_xlocal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xlyl.y, np.array(
                [-6.6742992, -5.8400118, -5.0057244, -4.171437, -3.3371496, -2.5028622, -1.6685748,
                 -0.8342874, 0., 0.8342874, 1.6685748, 2.5028622, 3.3371496, 4.171437, 5.0057244,
                 5.8400118]))

    def test_kxgrid_xglobal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xgyl.kx_fftorder, np.array(
                [0., 0.0392699, 0.0785398, 0.1178097, 0.1570796, 0.1963495, 0.2356194, 0.2748894,
                 0.3141593, 0.3534292, 0.3926991, 0.431969, -0.4712389, -0.43196899, -0.39269908,
                 -0.35342917, -0.31415927, -0.27488936, -0.23561945, -0.19634954, -0.15707963,
                 -0.11780972, -0.07853982, -0.03926991]))
        np.testing.assert_almost_equal(self.grid_xgyl.kx, np.array(
                [-0.4712389, -0.431969, -0.3926991, -0.3534292, -0.3141593, -0.2748894, -0.2356194,
                 -0.1963495, -0.1570796, -0.1178097, -0.0785398, -0.0392699, 0., 0.0392699,
                 0.0785398, 0.1178097, 0.1570796, 0.1963495, 0.2356194, 0.2748894, 0.3141593,
                 0.3534292, 0.3926991, 0.431969]))

    def test_kxgrid_xglobal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xgyg.kx_fftorder, np.array(
                [0., 0.07853982, 0.15707963, 0.23561945, 0.31415927, 0.39269908, 0.4712389,
                 0.54977871, 0.62831853, 0.70685835, 0.78539816, 0.86393798, -0.9424778,
                 -0.86393798, -0.78539816, -0.70685835, -0.62831853, -0.54977871, -0.4712389,
                 -0.39269908, -0.31415927, -0.23561945, -0.15707963, -0.07853982]))
        np.testing.assert_almost_equal(self.grid_xgyg.kx, np.array(
                [-0.9424778, -0.86393798, -0.78539816, -0.70685835, -0.62831853, -0.54977871,
                 -0.4712389, -0.39269908, -0.31415927, -0.23561945, -0.15707963, -0.07853982, 0.,
                 0.07853982, 0.15707963, 0.23561945, 0.31415927, 0.39269908, 0.4712389, 0.54977871,
                 0.62831853, 0.70685835, 0.78539816, 0.86393798]))

    def test_kxgrid_xlocal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xlyg.kx, np.array(
                [0., 0.03926991, 0.07853982, 0.11780972, 0.15707963, 0.19634954, 0.23561945,
                 0.27488936, 0.31415927, 0.35342917, 0.39269908, 0.43196899, 0.4712389, 0.51050881,
                 0.54977871, 0.58904862, 0.62831853, 0.66758844, 0.70685835, 0.74612826, 0.78539816,
                 0.82466807, 0.86393798, 0.90320789]))

    def test_kxgrid_xlocal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xlyl.kx_fftorder, np.array(
                [0., 0.03926991, 0.07853982, 0.11780972, 0.15707963, 0.19634954, 0.23561945,
                 0.27488936, 0.31415927, 0.35342917, 0.39269908, 0.43196899, 0.4712389, -0.43196899,
                 -0.39269908, -0.35342917, -0.31415927, -0.27488936, -0.23561945, -0.19634954,
                 -0.15707963, -0.11780972, -0.07853982, -0.03926991]))
        np.testing.assert_almost_equal(self.grid_xlyl.kx, np.array(
                [-0.43196899, -0.39269908, -0.35342917, -0.31415927, -0.27488936, -0.23561945,
                 -0.19634954, -0.15707963, -0.11780972, -0.07853982, -0.03926991, 0., 0.03926991,
                 0.07853982, 0.11780972, 0.15707963, 0.19634954, 0.23561945, 0.27488936, 0.31415927,
                 0.35342917, 0.39269908, 0.43196899]))

    def test_kygrid_xglobal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xgyl.ky, np.array(
                [0., 0.4707, 0.9414, 1.4121, 1.8828, 2.3535, 2.8242, 3.2949]))

    def test_kygrid_xglobal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xgyg.ky, np.array(
                [-1.8828, -1.4121, -0.9414, -0.4707, 0., 0.4707, 0.9414, 1.4121]))

    def test_kygrid_xlocal_yglobal(self):
        np.testing.assert_almost_equal(self.grid_xlyg.ky, np.array(
                [-1.8828, -1.4121, -0.9414, -0.4707, 0., 0.4707, 0.9414, 1.4121]))

    def test_kygrid_xlocal_ylocal(self):
        np.testing.assert_almost_equal(self.grid_xlyl.ky, np.array(
                [0., 0.4707, 0.9414, 1.4121, 1.8828, 2.3535, 2.8242, 3.2949]))

    def test_zgrid(self):
        np.testing.assert_almost_equal(self.grid_xgyl.z, np.array(
                [-3.1415927, -2.7488936, -2.3561945, -1.9634954, -1.5707963, -1.1780972, -0.7853982,
                 -0.3926991, 0., 0.3926991, 0.7853982, 1.1780972, 1.5707963, 1.9634954, 2.3561945,
                 2.7488936]))


class TestDiagSpace(unittest.TestCase):

    def setUp(self):
        cm = CommonData(".dat", -1, -2)
        self.grid = cm.spatialgrid

    def test_natural_fullspace(self):
        diagspace_natural = DiagSpace(self.grid, False, True, (0, None), (0, None), (0, None))
        np.testing.assert_array_equal(self.grid.x[diagspace_natural.xslice], self.grid.x[slice(24)])
        np.testing.assert_array_equal(self.grid.ky[diagspace_natural.yslice], self.grid.ky[slice(8)])
        np.testing.assert_array_equal(self.grid.z[diagspace_natural.zslice], self.grid.z[slice(16)])

    def test_natural_fullspace2(self):
        diagspace_natural = DiagSpace(self.grid, False, True, (None, None), (None, None), (None, None))
        np.testing.assert_array_equal(self.grid.x[diagspace_natural.xslice], self.grid.x[slice(24)])
        np.testing.assert_array_equal(self.grid.ky[diagspace_natural.yslice], self.grid.ky[slice(8)])
        np.testing.assert_array_equal(self.grid.z[diagspace_natural.zslice], self.grid.z[slice(16)])

    def test_natural_limitedspace(self):
        diagspace_natural = DiagSpace(self.grid, False, True, (5, -4), (2, -4), (5, -4))
        np.testing.assert_array_equal(self.grid.x[diagspace_natural.xslice], self.grid.x[slice(5, -4)])
        np.testing.assert_array_equal(self.grid.ky[diagspace_natural.yslice], self.grid.ky[slice(2, -4)])
        np.testing.assert_array_equal(self.grid.z[diagspace_natural.zslice], self.grid.z[slice(5, -4)])

    def test_natural_invalidspace(self):
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, False, True, (25, -4), (2, -4), (5, -4))
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, False, True, (1, -4), (-1, -4), (5, -4))
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, False, True, (1, -4), (2, -4), (4, 4))

    def test_xfourier_fullspace(self):
        diagspace = DiagSpace(self.grid, True, True, (0, None), (0, None), (0, None))
        np.testing.assert_array_equal(self.grid.kx[diagspace.xslice], self.grid.kx[slice(24)])
        np.testing.assert_array_equal(self.grid.ky[diagspace.yslice], self.grid.ky[slice(8)])
        np.testing.assert_array_equal(self.grid.z[diagspace.zslice], self.grid.z[slice(16)])

    def test_xfourier_limitedspace(self):
        diagspace = DiagSpace(self.grid, True, True, (5, -4), (2, -4), (5, -4))
        np.testing.assert_array_equal(self.grid.kx[diagspace.xslice], self.grid.kx[slice(5, -4)])
        np.testing.assert_array_equal(self.grid.ky[diagspace.yslice], self.grid.ky[slice(2, -4)])
        np.testing.assert_array_equal(self.grid.z[diagspace.zslice], self.grid.z[slice(5, -4)])

    def test_xfourier_invalidspace(self):
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, True, True, (25, -4), (2, -4), (5, -4))
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, True, True, (1, -4), (-1, -4), (5, -4))
        with np.testing.assert_raises(IndexError):
            DiagSpace(self.grid, True, True, (1, -4), (2, -4), (4, 4))