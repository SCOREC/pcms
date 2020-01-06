import unittest
import numpy as np

from pydiag.utils.geom import Geometry
from pydiag.utils.ParIO import Parameters
from pydiag.utils.comm import CommonData


class TestGeometry(unittest.TestCase):
    def setUp(self):
        par = Parameters()
        par.Read_Pars("parameters.dat")
        cm = CommonData(".dat", -1, -2)
        self.geom = Geometry(cm)

    def test_gxx(self):
        np.testing.assert_almost_equal(self.geom.gxx, np.ones((16, 24)))

    def test_dpdx_pm_arr(self):
        np.testing.assert_almost_equal(self.geom.dpdx_pm_arr, np.zeros(24))

    def test_q_prof(self):
        np.testing.assert_almost_equal(self.geom.q, np.array(
                [0.8559488, 0.8682008, 0.88713028, 0.91273725, 0.9450217, 0.98398364, 1.02962307,
                 1.08193999, 1.14093439, 1.20660628, 1.27895565, 1.35798251, 1.44368686, 1.53606869,
                 1.63512802, 1.74086482, 1.85327912, 1.9723709, 2.09814016, 2.23058692, 2.36971116,
                 2.51551289, 2.6679921, 2.8271488]))

    def test_jaco3d(self):
        np.testing.assert_almost_equal(np.sum(self.geom.jaco3d), 4303.3641548671203)


class TestGeometry3d(unittest.TestCase):

    def test_gxx(self):
        cm = CommonData("_3d", -1, -2)
        self.assertRaises(NotImplementedError, Geometry, cm)
