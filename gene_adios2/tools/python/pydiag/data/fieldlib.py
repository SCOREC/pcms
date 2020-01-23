from os.path import getsize

import numpy as np
from pydiag.data.base_file import BinaryFile, TimeSeries
from pydiag.utils.averages import z_av3d
from pydiag.utils.geom import Geometry


class FieldFile(BinaryFile):
    """ Class to parse binary field.dat diagnostic output of GENE

    :param filename: Should be "field.dat" or "field_fileextension"
    :param common: CommonData object of the run to analyse
    """
    def __init__(self, filename, common):
        super().__init__(filename, common)
        self._set_datatypes()
        self._set_sizes()
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
        self.file.seek(self.offset(var))
        var3d = np.fromfile(self.file, count=self.cm.pnt.nx0*self.cm.pnt.nky0*self.cm.pnt.nz0,
                            dtype=self.npct)
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


class FieldAverage(TimeSeries):
    """ Create an flux-surface averaged time series from the field file

    :param common: CommonData object of the run
    """

    def __init__(self, common):
        super().__init__("fieldaverage", common)
        self.cm = common
        self.phi0_ext = 0  # For the sinusoidal external potential
        self.nk_ext = 1
        self.k_ext = 1
        self.phase_phi_ext = 0
        self.Cxy = 1  # Geometry correction

    def generate_timeseries(self):
        """ Wrapper for fieldlib: Read the ky=0 component of the potential"""
        if not self.cm.y_local:
            raise NotImplementedError('y global is not supported!')
        try:
            self.phi0_ext = self.cm.pnt.phi0_ext
            self.nk_ext = self.cm.pnt.kxind_phi_ext
            self.phase_phi_ext = self.cm.pnt.phase_phi_ext
            self.k_ext = self.nk_ext/self.cm.pnt.lx*2*np.pi
        except AttributeError:
            self.phi0_ext = 0
            self.nk_ext = 1
            self.phase_phi_ext = 0
            self.k_ext = 1
        phifile = FieldFile("field" + self.cm.fileextension, self.cm)
        self.fileobject = phifile
        self.check_times()
        pos = self.calc_positions()
        self.timearray = np.array(self.timearray)[pos]
        geom = Geometry(self.cm)
        self.Cxy = geom.Cxy
        phi_ky0_read = np.zeros((len(self.timearray), self.cm.pnt.nx0), dtype=np.complex128)
        for tind, time in enumerate(self.timearray):
            phifile.set_time(time)
            phi_ky0_read[tind, :] = z_av3d(phifile.phi(), geom)[:, 0]
        self.dataarray = phi_ky0_read


def gluefieldaverage(fieldlist):
    """ Function to concatenate a list of field data into one single object

    for continuation runs
    :param fieldlist: List of FieldData objects to connect
    :returns: One connected FieldData object
    """
    # TODO: Consistency checks if runs match
    if len(fieldlist) == 1:
        return fieldlist[0]
    # Use first FieldFile as basis for the construction of the glued object
    result = fieldlist.pop(0)
    for field in fieldlist:
        # When times overlap, give following FieldFile preference
        while result.timearray[-1] > field.timearray[0]:
            del result.timearray[-1]
            del result.dataarray[-1]
        result.timearray = np.concatenate((result.timearray, field.timearray))
        result.phi = np.concatenate((result.dataarray, field.dataarray))
    result.endtime = fieldlist[-1].endtime
    del fieldlist
    return result
