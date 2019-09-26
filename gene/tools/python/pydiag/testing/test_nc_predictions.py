import unittest

import numpy as np

from pydiag.data.profile_data import ProfileData
from pydiag.utils import nc_predictions
from pydiag.utils.comm import CommonData
from pydiag.data.datafiles import RunDataFiles


class TestNCFuncs(unittest.TestCase):
    def setUp(self):
        """ Read a test run"""
        cm = CommonData(".dat", -2, -2)
        self.prof = ProfileData(cm, RunDataFiles(cm))
        self.prof.get_profiles()

    def test_changstandard(self):
        qch = self.prof.Qcharr
        chref = np.array([[ 0.263578,  0.160584,  0.127331,  0.11376,  0.108291,  0.106603,
         0.106457,  0.106367,  0.105245,  0.102356,  0.097369,  0.09039,
         0.081878,  0.072500,  0.062938,  0.053760,  0.045345,  0.037889,
         0.031447,  0.025976,  0.021385,  0.017565,  0.014391,  0.011801]])
        np.testing.assert_array_almost_equal(qch, chref)

    def test_calc_potatowidth(self):
        act_potwidth = nc_predictions.calc_potwidth(self.prof)
        ref_potwidth = 0.036087
        np.testing.assert_array_almost_equal(act_potwidth, ref_potwidth)

    def test_fb_param_hinton_haz(self):
        act_fbparam = nc_predictions.fb_param_hinton_haz(self.prof)
        ref_fbparam = np.array(
                [1.125644, 1.140127, 1.146389, 1.149796, 1.151817, 1.153022, 1.153675, 1.153913,
                 1.153811, 1.153416, 1.152762, 1.151883, 1.150819, 1.149616, 1.148326, 1.147004,
                 1.145699, 1.144451, 1.143285, 1.142214, 1.14124, 1.140356, 1.139548, 1.138801])
        np.testing.assert_array_almost_equal(act_fbparam, ref_fbparam)

    def test_fb_param_hinton_hirsch_sig(self):
        act_fbparam = nc_predictions.fb_param_hirsch_sig(self.prof)
        ref_fbparam = np.array(
                [0.967218, 0.917863, 0.864731, 0.811853, 0.759479, 0.707363, 0.655205, 0.602735,
                 0.549713, 0.49593, 0.441194, 0.385322, 0.328135, 0.269449, 0.209064, 0.146768,
                 0.082326, 0.015488, -0.054016, -0.126474, -0.202189, -0.28149, -0.364727,
                 -0.452287])
        np.testing.assert_array_almost_equal(act_fbparam, ref_fbparam)


if __name__ == '__main__':
    unittest.main()
