""" Contains the class to read binary srcmom files"""

from os.path import getsize
import numpy as np

from pydiag.data.base_file import BinaryFile, TimeSeries


class SrcmomFile(BinaryFile):
    """ Class to read a binary srcmom output file

    A srcmom file entry consists of the time stamp followed by radial profiles for (1)particle,
     (2)momentum and (3)heat input by (1)the Krook heat source, (2)the Krook particle source and (3)
     the f0 contributions (both the localized source and the neoclassical term)

    :param filename: Typically srcmom_spec_fext
    :param common: CommonData object of the run
    """

    def __init__(self, filename, common, usemmap=False):
        super().__init__(filename, common, usemmap=usemmap)
        self.prev_data = (0, "", "")
        self.filename = filename
        self.srcdata = None
        self._set_sizes()

    def _set_sizes(self):
        """ Set up the sizes for the records in Fortran binary files"""
        super()._set_sizes()
        self.nmoms = 3
        self.nsrc = 3
        self.moms = {"part": 0, "mom": 1, "energy": 2}
        self.sources = {"ck_heat": 0, "ck_part": 1, "f0": 2}
        # Each dataset at a time consists of 3 records (the source types) which contain 3 nx arrays
        # (particle, heat, momentum moment). Each record is surrounded by its length as a 4 byte int
        self.entrysize = self.cm.pnt.nx0*self.realsize
        self.leapsrcmom = self.nmoms*(self.nsrc*self.entrysize + 2*self.intsize)

    def get_timearray(self):
        """ Get the time array from the file"""
        # Fortran has a particular binary format: see time_entry()
        for i in range(int(getsize(self.filename)/(self.leapsrcmom + self.tesize))):
            self.timearray.append(float(self.te.unpack(self.file.read(self.tesize))[1]))
            self.file.seek(self.leapsrcmom, 1)

    def reset_tinds(self):
        """ Reset time index and previously read data """
        super().reset_tinds()
        self.prev_data = (0, "", "")

    def readvar(self, var, subvar):
        """ Read one of the nine vectors at a given time

        :param var: one of self.moms values
        :param subvar: one of self.sources values
        """
        if not (self.tind, var, subvar) == self.prev_data:
            self.srcdata = self.readfunc(self.file, count=self.cm.pnt.nx0, dtype=self.nprt,
                                         offset=self.offset(var, subvar))
            self.prev_data = (self.tind, var, subvar)
        return self.srcdata

    def offset(self, var, subvar):
        """ Calculate offset in field file for a given timestep and variable"""
        return self.tind*(self.tesize + self.leapsrcmom) + self.tesize + self.intsize + subvar*(
            self.nmoms*self.entrysize + 2*self.intsize) + var*self.entrysize

    def srcmoment(self, momtype, srctype):
        """ Return the source moment, read from the file if necessary

        :param momtype: Type of the moment (key of self.moms)
        :param srctype: Type of the source (key of self.source)
        """
        try:
            imom = self.moms[momtype]
            isrc = self.sources[srctype]
        except KeyError:
            raise KeyError("Not a valid source or source moment entry")
        return self.readvar(imom, isrc)

    def get_var(self, varname):
        mom, src = varname.split("_")
        return self.srcmoment(mom, src)

class SrcmomSeries(TimeSeries):
    """ Construct a time series of the source moment diag """

    def __init__(self, common, rundatafiles, momtype, srctype):
        super().__init__("srcmomseries", common)
        self.momtype = momtype
        self.srctype = srctype
        self.srcmoms = [[] for _ in common.specnames]
        self.rundatafiles = rundatafiles

    def generate_timeseries(self):
        srcmomfiles = [self.rundatafiles.get_fileobject("srcmom_{}".format(spec))
                       for spec in self.cm.specnames]
        self.fileobject = srcmomfiles[0]
        self.check_times()
        tinds = self.calc_positions()
        for ispec, smf in enumerate(srcmomfiles):
            for pos in tinds:
                smf.set_tind(pos)
                self.srcmoms[ispec].append(np.array(smf.srcmoment(self.momtype, self.srctype)))
        self.timearray = np.array(self.timearray)[tinds]