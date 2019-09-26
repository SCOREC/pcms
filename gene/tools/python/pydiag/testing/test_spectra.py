import unittest
import numpy as np

from pydiag.data.datafiles import RunDataFiles
from pydiag.data.slices import MomFieldSlice
from pydiag.utils.comm import CommonData, DiagSpace
from pydiag.diagplots.plot_spectra import AmplitudeSpectra, FluxSpectra


class TestAmplitudeSpectra(unittest.TestCase):

    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.rds = RunDataFiles(self.cm)

    def test_phi_spec_ky(self):
        ampspec = AmplitudeSpectra(self.cm, "ions", self.rds, moms=["phi"])
        np.testing.assert_almost_equal(ampspec.momamps_ky["phi"].timeaverage, np.array(
                [7.4345451e-09, 5.1101912e-03, 5.2098234e-04, 2.3990140e-05, 5.0293696e-07,
                 5.1400817e-09, 2.8572088e-11, 9.7016507e-14]))

    def test_tpar_spec_ky(self):
        ampspec = AmplitudeSpectra(self.cm, "ions", self.rds, moms=["tpar"])
        np.testing.assert_almost_equal(ampspec.momamps_ky["tpar"].timeaverage, np.array(
                [2.7386640e-10, 2.0322812e-04, 7.7205893e-06, 5.7506578e-06, 3.3685789e-07,
                 6.7326446e-09, 6.2655218e-11, 3.1763233e-13]))

    def test_phi_spec_kx(self):
        ampspec = AmplitudeSpectra(self.cm, "ions", self.rds, moms=["phi"])
        print(repr(ampspec.momamps_kx["phi"].timeaverage))
        np.testing.assert_almost_equal(ampspec.momamps_kx["phi"].timeaverage, np.array(
                [2.32447438e+01, 1.27613993e+01, 2.44130432e+00, 1.57438019e-01, 3.38981831e-03,
                 2.89089079e-05, 2.25615905e-07, 4.48767729e-09, 1.96479193e-08, 6.92885159e-08,
                 5.39650901e-08, 2.05277960e-08, 4.56397473e-09, 2.51264624e-08, 4.90561181e-08,
                 3.91471289e-08, 2.50360322e-08, 2.24943989e-09, 8.53812082e-07, 1.23777063e-04,
                 9.40146514e-03, 2.98067240e-01, 3.52417607e+00, 1.50896683e+01]))
