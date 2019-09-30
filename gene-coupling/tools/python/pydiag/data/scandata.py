# -*- coding: utf-8 -*-
""" Module to handle the output of a scan

including local scans derived from global profiles and scans producing neoclassical data
"""
import numpy as np

from pydiag.utils.ParIO import Parameters


class Scandata(object):
    """ Collect the data of a scan and store them in nice numpy arrays"""
    def __init__(self):

        mainpars = Parameters()
        mainpars.Read_Pars("parameters")
        try:
            # The dimensions in the namelist are in inverse order of the scan.log file!
            self.scandims = [int(dim) for dim in mainpars.pardict["scan_dims"].split()][::-1]
        except AttributeError:
            # This occurs for only one scan dimension
            self.scandims = [mainpars.pardict["scan_dims"]]
        # Build a field for each scan parameter with the full scandimensions (concatenating tuples)
        self.grid = np.empty((tuple(self.scandims)+(len(self.scandims),)))
        # its also nice to have them as a list of 1d vectors
        self.redgrid = [np.empty(dim) for dim in self.scandims]
        self.growthrates = np.empty(self.scandims)
        self.frequencies = np.empty(self.scandims)
        self.neodata = {"Gammanc": np.empty(self.scandims), "Qnc": np.empty(self.scandims),
                        "Pinc": np.empty(self.scandims), "jbs": np.empty(self.scandims)}
        self.scannames = []     # Names of the parameters scanned
        self.fetchedpars = {}   # Contains the data fetched from individual parameter_xxxx files

    def parsescan(self):
        """ Parse scan.log and get additional data from parameters_run """
        try:
            with open("scan.log", 'r') as scanlog:
                header = scanlog.readline().split()
                self.set_scan_names(header)
                for line in scanlog:
                    lst = line.split()
                    # scan.log has column majority ("Fortran") order
                    coord = np.unravel_index(int(lst[0])-1, self.scandims, order='F')
                    for ipar in range(len(self.scandims)):
                        self.grid[coord + (ipar,)] = float(lst[2+2*ipar])
                    self.growthrates[coord] = float(lst[-2])
                    self.frequencies[coord] = float(lst[-1])
        except IOError:
            print("Could not read scan.log file")
            raise
        self._reducegrid()

    def _reducegrid(self):
        """ Reduce the big grid variable for the scan parameters to a list of arrays """
        if len(self.scandims) == 1:
            self.redgrid = {self.scannames: self.grid}
        else:
            for idim, dimsize in enumerate(self.scandims):
                self.redgrid[idim] = np.rollaxis(self.grid[..., idim], idim).flatten(order='F')[
                                     0:dimsize]
            self.redgrid = {self.scannames[i]: self.redgrid[i] for i, _ in enumerate(self.scannames)}

    def parseneo(self):
        """ Parse neo.log from simulations with neoclassical diagnostics """
        self.parsescan()
        try:
            with open("neo.log", 'r') as neolog:
                neolog.readline()
                for line in neolog:
                    lst = line.split()
                    # neo.log has column majority ("Fortran") order
                    coord = self._scannum_to_coord(int(lst[0]))
                    self.neodata["Gammanc"][coord] = float(lst[-4])
                    self.neodata["Qnc"][coord] = float(lst[-3])
                    self.neodata["Pinc"][coord] = float(lst[-2])
                    self.neodata["jbs"][coord] = float(lst[-1])
        except IOError:
            print("Could not read neo.log file. Is this a neoclassical run?")
            raise

    def _scannum_to_coord(self, scannum):
        return np.unravel_index(scannum-1, self.scandims, order='F')

    def fetchpars(self, pars, reread=False):
        """ Get parameter values from the parameters_fext files

        Warning: This can take quite long depending on the file
        system and scan dimensions
        :param pars: List of the parameters to read
        :param reread: Force update of the already parsed parameters
        """
        if not pars:
            return
        if (not reread) and (self.fetchedpars.keys() == pars):
            return
        runpar = Parameters()
        self.fetchedpars = {parname: np.full(self.scandims, np.nan) for parname in pars}
        # The files are numbered from one to scandim.size
        for runnum in range(1, self.growthrates.size+1):
            fext = str(runnum).zfill(4)
            runpar.Read_Pars("parameters_"+fext)
            coord = self._scannum_to_coord(runnum)
            par_not_found = []
            for parkey in pars:
                try:
                    self.fetchedpars[parkey][coord] = runpar.pardict[parkey]
                except KeyError:
                    par_not_found.append(parkey)
            if par_not_found:
                print("Parameter(s) {} not found. Will return NaN for it.".format(par_not_found))
                for rp in par_not_found:
                    pars.remove(rp)
                par_not_found.clear()

    def set_scan_names(self, header):
        """ Extract scan parameter names

        :param header: First line of scan.log or neo.log files
        """
        self.scannames.clear()
        for ipar in range(len(self.scandims)):
            # First comes "#Run" then |, then the parameter, a spec number
            # another | until / followed by /Eigenvalue1
            self.scannames.append(header[2 + 3*ipar])
