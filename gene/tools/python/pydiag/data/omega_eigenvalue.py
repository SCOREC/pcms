# -*- coding: utf-8 -*-
""" Module for reading omega and eigenvalue output files """
import numpy as np

class Eigenvaluedata(object):
    """ Collect the data of an oemga or eigenvalue file and store it in a nice numpy array"""

    def __init__(self,common):
        if common.pnt.nonlinear:
            raise RuntimeError("Not a linear simulation - won't read eigenvalues")

        if common.pars["comp_type"]=="'IV'":
            self.parseomegafile(common.fileextension)
        else:
            self.parseeigenvaluesfile(common.fileextension)

    def parseomegafile(self,fileextension):
        """ Parse omega.dat """
        try:
            with open("omega"+fileextension, 'r') as omegafile:
                data = omegafile.readline().split()
                self.growth_rate = float(data[-2])
                self.frequency = float(data[-1])
        except IOError:
            print("Could not read omega"+fileextension+" file")
            raise

    def parseeigenvaluesfile(self):
        """ Parse eigenvalues.dat """
        try:
            with open("eigenvalues"+fileextension, 'r') as evfile:
                data = evfile.readline().split()
                self.growth_rate = float(data[-2])
                self.frequency = float(data[-1])
        except IOError:
            print("Could not read eigenvalues"+fileextension+" file")
            raise
        raise NotImplementedError
