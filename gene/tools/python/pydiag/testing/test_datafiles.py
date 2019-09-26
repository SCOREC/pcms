import unittest
import numpy as np

from pydiag.utils.comm import CommonData, SpatialGrid, DiagSpace
from pydiag.data.datafiles import RunDataFiles


class TestRunDataFiles(unittest.TestCase):

    def setUp(self):
        self.cm = CommonData(".dat", -1, -2)
        self.rdf = RunDataFiles(self.cm)

    def test_get_fileobject_mom(self):
        mom = self.rdf.get_fileobject("mom_ions")
        np.testing.assert_almost_equal(np.sum(mom.tpar()),
                                       -0.046149018431413114+1.6268003468905556e-05j)

    def test_get_fileobject_field(self):
        field = self.rdf.get_fileobject("field")
        np.testing.assert_almost_equal(np.sum(field.phi()),
                                       (29.853139580890563 - 7.0717609891371436e-06j))
        field2 = self.rdf.get_fileobject("field")
        np.testing.assert_almost_equal(np.sum(field2.phi()),
                                       (29.853139580890563 - 7.0717609891371436e-06j))


    def test_get_fileobject_nrg(self):
        nrg = self.rdf.get_fileobject("nrg")
        nrg.generate_timeseries()
        np.testing.assert_almost_equal(np.sum(nrg.dataarray), 0.5308799971974911)

    def test_get_fileobject_error(self):
        with np.testing.assert_raises(RuntimeError):
            self.rdf.get_fileobject("foobar")
