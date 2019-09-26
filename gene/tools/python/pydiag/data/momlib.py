import numpy as np
from os.path import getsize

from pydiag.data.base_file import BinaryFile


class MomFile(BinaryFile):
    """ Class to parse binary mom_[spec].dat diagnostic output of GENE

    :param filename: mom file name
    :param common: CommonData object of the run to analyse
    """
    def __init__(self, filename, common, usemmap=False):
        super().__init__(filename, common, usemmap=usemmap)
        self.momvar3d = {"dens": None, "tpar": None, "tperp": None, "qpar": None, "qperp": None,
                         "upar": None, "n00": None, "n20": None, "n02": None}

    def _set_sizes(self):
        super()._set_sizes()
        self.nmoms = self.cm.pnt.n_moms
        self.entrysize = self.cm.pnt.nx0*self.cm.pnt.nky0*self.cm.pnt.nz0*self.complexsize
        # jumps in bytes in mom files
        self.leapmom = self.nmoms*(self.entrysize + 2*self.intsize)

    def reset_tinds(self):
        """ Reset time indices when reloading file"""
        super().reset_tinds()
        self.mtind = [-1]*self.nmoms

    def get_timearray(self):
        """Get time array for mom file """
        self.timearray = []
        self.file.seek(0)
        for i in range(int(getsize(self.filename)/(self.leapmom + self.tesize))):
            self.timearray.append(float(self.te.unpack(self.file.read(self.tesize))[1]))
            self.file.seek(self.leapmom, 1)

    def offset(self, var):
        """Calculate offset in mom file for a given self.time and variable"""
        if var in range(self.nmoms):
            leap = self.leapmom
            return self.tesize + self.tind*(self.tesize + leap) + var*(
                self.entrysize + 2*self.intsize) + self.intsize

    def readvar(self, var):
        """ Return 3d mom data at the time set in self.time"""
        var3d = self.readfunc(self.file, count=self.cm.pnt.nx0*self.cm.pnt.nky0*self.cm.pnt.nz0,
                            dtype=self.npct, offset=self.offset(var))
        if self.cm.x_local and not self.cm.y_local:  # y-global has yx order
            var3d = var3d.reshape(self.cm.pnt.nky0, self.cm.pnt.nx0, self.cm.pnt.nz0, order="F")
            var3d = np.swapaxes(var3d, 0, 1)
        else:
            var3d = var3d.reshape(self.cm.pnt.nx0, self.cm.pnt.nky0, self.cm.pnt.nz0, order="F")
        return var3d

    def get_var(self, varname):
        """ Return a 3d moment by name"""
        varidx = {"dens": 0, "tpar": 1, "tperp": 2, "qpar": 3, "qperp": 4, "upar": 5, "n00": 6, "n20": 7, "n02": 8}
        if not self.mtind[varidx[varname]] == self.tind:
                self.mtind[varidx[varname]] = self.tind
                self.momvar3d[varname] = self.readvar(varidx[varname])
        return self.momvar3d[varname]

    def dens(self):
        """ Return density fluctuation"""
        if not self.mtind[0] == self.tind:
            self.mtind[0] = self.tind
            self.dens3d = self.readvar(0)
        return self.dens3d

    def tpar(self):
        """ Return T_\parallel fluctuations"""
        if not self.mtind[1] == self.tind:
            self.mtind[1] = self.tind
            self.tpar3d = self.readvar(1)
        return self.tpar3d

    def tperp(self):
        """ Return T_\perp fluctuations"""
        if not self.mtind[2] == self.tind:
            self.mtind[2] = self.tind
            self.tperp3d = self.readvar(2)
        return self.tperp3d

    def qpar(self):
        """ Return q_\parallel fluctuations"""
        if not self.mtind[3] == self.tind:
            self.mtind[3] = self.tind
            self.qpar3d = self.readvar(3)
        return self.qpar3d

    def qperp(self):
        """ Return q_\perp fluctuations"""
        if not self.mtind[4] == self.tind:
            self.mtind[4] = self.tind
            self.qperp3d = self.readvar(4)
        return self.qperp3d

    def upar(self):
        """ Return parallel flow fluctuations"""
        if not self.mtind[5] == self.tind:
            self.mtind[5] = self.tind
            self.upar3d = self.readvar(5)
        return self.upar3d

    def n00(self):
        """ Return perpendicular pressure fluctuations"""
        if not self.mtind[6] == self.tind:
            self.mtind[6] = self.tind
            self.n00_3d = self.readvar(6)
        return self.n00_3d

    def n20(self):
        """ Return N20 fluctuations"""
        if not self.mtind[7] == self.tind:
            self.mtind[7] = self.tind
            self.n20_3d = self.readvar(7)
        return self.n20_3d

    def n02(self):
        """ Return N02 fluctuations"""
        if not self.mtind[8] == self.tind:
            self.mtind[8] = self.tind
            self.n02_3d = self.readvar(8)
        return self.n02_3d
