import unittest
import numpy as np
import os

from pydiag.utils.comm import CommonData
from pydiag.data import profile_data as prd


# TODO: Test for initial profiles, nustar
class TestProfileData(unittest.TestCase):
    def setUp(self):
        """ Read a test run"""
        self.cm = CommonData(".dat", -1, -2)
        self.prof = prd.ProfileData(self.cm)

    def test_timefield(self):
        """ Make sure the time stamps are read correctly"""
        self.prof.get_profiles()
        times = self.prof.prspec[0].timearray
        np.testing.assert_almost_equal(times, np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538]))

    def test_qprofile(self):
        """ Check if the safety factor is correctly read

        Depends on the geometry file parser
        """
        np.testing.assert_almost_equal(self.prof.q,
                                       [0.8559488, 0.8682008, 0.88713028, 0.91273725, 0.9450217,
                                        0.98398364, 1.02962307, 1.08193999, 1.14093439, 1.20660628,
                                        1.27895565, 1.35798251, 1.44368686, 1.53606869, 1.63512802,
                                        1.74086482, 1.85327912, 1.9723709, 2.09814016, 2.23058692,
                                        2.36971116, 2.51551289, 2.6679921, 2.8271488])

    def test_calc_afs(self):
        self.prof.calc_afs()
        np.testing.assert_almost_equal(self.prof.Afs,
                                       [0.7390386, 1.2926996, 1.8463606, 2.4000216, 2.9536827,
                                        3.5073437, 4.0610047, 4.6146657, 5.1683267, 5.7219877,
                                        6.2756488, 6.8293098, 7.3829708, 7.9366318, 8.4902928,
                                        9.04395386, 9.59761487, 10.15127589, 10.70493691,
                                        11.25859793, 11.81225894, 12.36591996, 12.91958098,
                                        13.473242])

    def test_generatetimetrace(self):
        self.prof.get_profiles()
        self.prof.generate_timetrace(0.5)
        ttdata = np.loadtxt("ttrace_0.5Qturb_ions.dat")
        np.testing.assert_almost_equal(ttdata, np.array(
                [[0.00000000e+00, -1.90706605e-11], [1.07600000e-01, 2.06563391e-02],
                 [2.15200000e-01, 4.09225683e-02], [3.22800000e-01, 6.04909977e-02],
                 [4.30400000e-01, 7.91618971e-02], [5.38000000e-01, 9.68863567e-02]]))
        os.remove("ttrace_0.5Qturb_ions.dat")

    def test_gluefunction(self):
        cm2 = CommonData(".dat", -2, -2)
        prof2 = prd.ProfileData(cm2)
        prof2.get_profiles()
        self.prof.get_profiles()
        glue = prd.glueprofiledata([self.prof, prof2])
        np.testing.assert_almost_equal(glue.prspec[0].timearray,
                                       np.array([0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538]))


class TestProfileFileData(unittest.TestCase):
    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.pfd = prd.ProfileFile("profile_ions.dat", self.cm)

    def test_invalidfile(self):
        with np.testing.assert_raises(IOError):
            temp = prd.ProfileFile("foobar", self.cm)

    def test_timefield(self):
        times = self.pfd.timearray
        np.testing.assert_almost_equal(times, [0., 0.1076, 0.2152, 0.3228, 0.4304, 0.538])

    def test_getminmaxtime(self):
        np.testing.assert_almost_equal(self.pfd.get_minmaxtime(), [0, 0.538])

    def test_singletime(self):
        cm2 = CommonData(".dat", 0.3, 0.3)
        temppfd = prd.ProfileFile("profile_ions.dat", cm2)
        temppfd.generate_timeseries()
        # The class should round to the nearest existing timestamp
        np.testing.assert_almost_equal(temppfd.timearray, [0.3228])

    def test_changetolasttime(self):
        """ Set the class to read only the last time point, then reset"""
        self.pfd.starttime = -2
        self.pfd.endtime = -2
        self.pfd.generate_timeseries()
        np.testing.assert_almost_equal(self.pfd.timearray, [0.538])
        self.pfd.starttime = -1
        self.pfd.endtime = -2
        self.pfd.generate_timeseries()

    def test_temperature_profile(self):
        temp = self.pfd.dataarray["Ts"]
        np.testing.assert_almost_equal(temp, [np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376995, 1.265918, 1.155566, 1.049943, 0.9524323, 0.8653768, 0.7899404, 0.7262189,
                 0.6735136, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852]), np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376995, 1.265918, 1.155566, 1.049943, 0.9524323, 0.8653769, 0.7899405, 0.7262189,
                 0.6735136, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852]), np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376995, 1.265918, 1.155566, 1.049943, 0.9524325, 0.8653771, 0.7899406, 0.726219,
                 0.6735136, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852]), np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376995, 1.265918, 1.155565, 1.049943, 0.9524328, 0.8653775, 0.7899409, 0.7262191,
                 0.6735136, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852]), np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376994, 1.265917, 1.155565, 1.049943, 0.9524332, 0.865378, 0.7899412, 0.7262193,
                 0.6735137, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852]), np.array(
                [1.972833, 1.933577, 1.885395, 1.827106, 1.757821, 1.677191, 1.585674, 1.484751,
                 1.376994, 1.265916, 1.155565, 1.049943, 0.9524337, 0.8653786, 0.7899417, 0.7262195,
                 0.6735137, 0.6306467, 0.5962352, 0.5688862, 0.5473135, 0.5303927, 0.5171762,
                 0.5068852])])
