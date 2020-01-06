import numpy as np
import pydiag.utils.averages as avg
from pydiag.utils.geom import Geometry
from pydiag.data.base_file import TimeSeries
from pydiag.data.fieldlib import FieldFile
from pydiag.data.momlib import MomFile


def abssquare(var):
    return np.abs(var) ** 2


class MomFieldSlice(TimeSeries):
    """ Class to generate spatial slices of Mom or Field data"""

    def __init__(self, common, varname, species, diagspace, modifier_func=abssquare):
        """ Constructor for MomFieldSlice

        :param common: CommonData object of the run
        :param varname: Name of the variable to process
        :param species: Name of the species to process
        :param diagspace: Fourier or configuration space slice desired for x and y, are
        transforms necessary from the native format of the run?
        :param modifier_func: A function reference which should be applied to the 3d data first
        By default the squared average is calculated
        """
        super().__init__("slice", common)
        self.cm = common
        self.geom = Geometry(self.cm)
        self.diagspace = diagspace
        self.varname = varname
        if varname in ["phi", "apar", "bpar"]:
            self.fileobject = FieldFile("field{}".format(common.fileextension), common)
        elif varname in ["dens", "tpar", "tperp", "qpar", "qperp", "upar"]:
            self.fileobject = MomFile("mom_{}{}".format(species, common.fileextension), common)
        else:
            raise RuntimeError("Invalid variable name")
        self.modifier_func = modifier_func

    def generate_timeseries(self):
        self.check_times()
        pos = self.calc_positions()
        self.dataarray = []
        new_timearray = []
        for tind in pos:
            self.dataarray.append(self.generate_slice_attime(self.timearray[tind]))
            new_timearray.append(self.timearray[tind])
        self.timearray = new_timearray
        self.dataarray = np.array(self.dataarray)

    def generate_slice_attime(self, time):
        """ Create a slice of the data at a fixed time with the desired dimensions

        :param time: The exact time to use
        """
        self.fileobject.set_time(time)
        var = self.fileobject.get_var(varname=self.varname)
        var = self.modifier_func(var)
        averagedvar = avg.av3d_by_switch(self.diagspace.xavg, self.diagspace.yavg,
                                         self.diagspace.zavg)(var, self.geom)
        return averagedvar[self.diagspace.diagslice]

    def write_to_numpy(self, filename):
        self.generate_timeaverage()
        desc = """Slice diagnostic for {} (with {} applied), t={}, xavg,yavg,zavg={},{},{},"""
        """Range (x,y,z)={},{},{}""".format(self.varname, self.modifier_func, self.get_minmaxtime(),
                                            self.diagspace.xavg, self.diagspace.yavg,
                                            self.diagspace.zavg, self.diagspace.xslice,
                                            self.diagspace.yslice, self.diagspace.zslice)
        np.savez_compressed(filename, description=desc, timearray=self.timearray,
                            dataarray=self.dataarray, timeaverage=self.timeaverage)

    def write_to_ascii(self, filename, t_average=False):
        if np.squeeze(self.dataarray[0, ...]).ndim > 2:
            raise RuntimeError("Unable to write 4d data (time + x,y,z) as ASCII file")
        if t_average:
            self.generate_timeaverage()
            np.savetxt(filename, np.squeeze(self.timeaverage),
                       header="Slice diagnostic for {} (with {} applied), averaged t={}\n".format(
                               self.varname, self.modifier_func, self.get_minmaxtime()))
        else:
            with open(filename, "wb") as file:
                file.write(bytearray("""# Slice diagnostic for {} (with {} applied), t={},\n"""
                                     """# xavg,yavg,zavg={},{},{}, Range (x,y,z)={},{},{}\n""".format(
                        self.varname, self.modifier_func, self.get_minmaxtime(),
                        self.diagspace.xavg, self.diagspace.yavg, self.diagspace.zavg,
                        self.diagspace.xslice, self.diagspace.yslice, self.diagspace.zslice),
                        encoding="utf-8"))
                for tind, time in enumerate(self.timearray):
                    np.savetxt(file, np.atleast_1d(np.squeeze(self.dataarray[tind])), fmt='%.9e',
                               header="{}".format(self.timearray[tind]))
                    file.write(b"\n")
