from os.path import getsize

import numpy as np
import pydiag.data.base_file


class FieldFile(pydiag.data.base_file.BinaryFile):
    """ Class to parse binary field.dat diagnostic output of GENE

    :param filename: Should be "field.dat" or "field_fileextension"
    :param common: CommonData object of the run to analyse
    """
    def __init__(self, filename, common, usemmap=False):
        super().__init__(filename, common, usemmap=usemmap)
        self.ftind = [-1]*self.nfields
        self.fieldvar3d = {"phi": None, "apar": None, "bpar": None}
        self.phi3d = None
        self.apar3d = None
        self.bpar3d = None

    def _set_sizes(self):
        super()._set_sizes()
        self.nfields = int(self.cm.pars['n_fields'])
        self.entrysize = self.cm.pnt.nx0*self.cm.pnt.nky0*self.cm.pnt.nz0*self.complexsize
        # jump in bytes in field files
        self.leapfld = self.nfields*(self.entrysize + 2*self.intsize)

    def get_timearray(self):
        """Get time array for field file """
        self.timearray = []
        self.file.seek(0)
        for i in range(int(getsize(self.filename)/(self.leapfld + self.tesize))):
            self.timearray.append(float(self.te.unpack(self.file.read(self.tesize))[1]))
            self.file.seek(self.leapfld, 1)

    def offset(self, var):
        """Calculate offset in field file for a given self.time and variable"""
        if var in [i for i in range(self.nfields)]:
            return self.tesize + self.tind*(self.tesize + self.leapfld) + var*(
                self.entrysize + 2*self.intsize) + self.intsize

    def readvar(self, var):
        """ Return 3d field data at the time set in self.time"""
        var3d = self.readfunc(self.file, count=self.cm.pnt.nx0*self.cm.pnt.nky0*self.cm.pnt.nz0,
                              dtype=self.npct, offset=self.offset(var))

        # Bring array into x, y. z order
        if self.cm.x_local and not self.cm.y_local:  # y-global has yx order
            var3d = var3d.reshape(self.cm.pnt.nky0, self.cm.pnt.nx0, self.cm.pnt.nz0, order="F")
            var3d = np.swapaxes(var3d, 0, 1)
        else:
            var3d = var3d.reshape(self.cm.pnt.nx0, self.cm.pnt.nky0, self.cm.pnt.nz0, order="F")
        return var3d

    def reset_tinds(self):
        """ Reset time indices when reloading file"""
        super().reset_tinds()
        self.ftind = [-1]*self.nfields

    def get_var(self, varname):
        """ """
        varidx = {"phi": 0, "apar": 1, "bpar": 2}
        if varidx[varname] < self.nfields:
            if not self.ftind[varidx[varname]] == self.tind:
                self.ftind[varidx[varname]] = self.tind
                self.fieldvar3d[varname] = self.readvar(varidx[varname])
        return self.fieldvar3d[varname]

    def phi(self):
        """ Return phi (electrostatic potential)

        Read from file if necessary
        """
        if not self.ftind[0] == self.tind:
            self.ftind[0] = self.tind
            self.phi3d = self.readvar(0)
        return self.phi3d

    def apar(self):
        """ Return apar (parallel component of magnetic vector potential)

        Read from file if necessary
        """
        if self.nfields > 1:
            if not self.ftind[1] == self.tind:
                self.ftind[1] = self.tind
                self.apar3d = self.readvar(1)
            return self.apar3d

    def bpar(self):
        """ Return bpar (parallel component of magnetic field fluctuation)

        Read from file if necessary
        """
        if self.nfields > 2:
            if not self.ftind[2] == self.tind:
                self.ftind[2] = self.tind
                self.bpar3d = self.readvar(2)
            return self.bpar3d
