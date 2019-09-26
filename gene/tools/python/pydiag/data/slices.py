import numpy as np
import warnings
import pydiag.utils.averages as avg
from pydiag.utils.geom import Geometry
from pydiag.data.base_file import TimeSeries
import pydiag.utils.fourier as fourier


class MomFieldSlice(TimeSeries):
    """ Class to generate spatial slices of Mom or Field data"""

    def __init__(self, common, quantity, species, diagspace, rundatafiles,
                 modifier_func=lambda var: np.real(var*np.conj(var))):
        """ Constructor for MomFieldSlice

        :param common: CommonData object of the run
        :param quantity: Name of the variable to process
        :param species: Name of the species to process
        :param diagspace: Fourier or configuration space slice desired for x and y, are
        transforms necessary from the native format of the run?
        :param rundatafiles: A RunDataFiles object that provides the mom and field files necessary
        :param modifier_func: A function reference which should be applied to the 3d data first
        By default the squared absolute is calculated
        """
        super().__init__("slice", common)
        self.cm = common
        self.geom = Geometry(self.cm)
        self.diagspace = diagspace
        self.quantity = quantity
        if quantity in ["phi", "apar", "bpar"]:
            self.fileobject = rundatafiles.get_fileobject("field")
        elif quantity in ["dens", "tpar", "tperp", "qpar", "qperp", "upar",
                          "n00", "n20", "n02"]:
            self.fileobject = rundatafiles.get_fileobject("mom_{}".format(species))
        else:
            raise RuntimeError("Invalid variable name")
        self.modifier_func = modifier_func

    def generate_timeseries(self, sparsefactor=1):
        self.check_times()
        pos = self.calc_positions()
        self.dataarray = []
        new_timearray = []
        if sparsefactor != 1:
            pos = pos[::sparsefactor]
        for tind in pos:
            self.dataarray.append(self.generate_slice_attime(self.timearray[tind]))
            new_timearray.append(self.timearray[tind])
        self.timearray = np.array(new_timearray)
        self.dataarray = np.atleast_2d((np.array(self.dataarray)))

    def generate_slice_attime(self, time):
        """ Create a slice of the data at a fixed time with the desired dimensions

        :param time: The exact time to use
        """
        self.fileobject.set_time(time)
        var = self.fileobject.get_var(varname=self.quantity)
        averagedvar = avg.av3d_by_switch(self.diagspace.xavg, self.diagspace.yavg,
                                         self.diagspace.zavg)(var, self.geom)
        averagedvar = self._apply_fouriertransforms(averagedvar)
        averagedvar = self.modifier_func(averagedvar)
        return averagedvar[self.diagspace.diagslice]

    def _apply_fouriertransforms(self, var):
        if self.diagspace.xavg and self.diagspace.yavg:  # Nothing to do if both are averaged
            return var
        xaxis = -3  # By default the last three dimensions of var are x, y, z.
        yaxis = -2  # This is changed if averages are applied
        zaxis = -1
        if self.diagspace.yavg:
            xaxis += 1
            yaxis += 1
        if self.diagspace.zavg:
            xaxis += 1
            yaxis += 1

        if self.diagspace.x_fourier != self.cm.x_local and not self.diagspace.xavg:
            if self.cm.x_local:
                var = fourier.kx_to_x(var, self.cm.pnt.nx0, axis=xaxis)
            else:
                var = fourier.x_to_kx(var, self.cm.pnt.nx0, axis=xaxis)
        if self.diagspace.y_fourier != self.cm.y_local and not self.diagspace.yavg:
            if self.cm.y_local:
                var = fourier.ky_to_y(var, self.cm.pnt.nky0, axis=yaxis)
            else:
                warnings.warn("y-global is not well tested", RuntimeWarning)
                var = fourier.y_to_ky(var, self.cm.pnt.nky0, axis=yaxis)
        if self.diagspace.z_fourier:
            var=fourier.z_to_kz(var, self.cm.pnt.nz0, axis=zaxis)
        return var

    def write_to_numpy(self, filename):
        self.generate_timeaverage()
        desc = """Slice diagnostic for {} (with {} applied), t={}, xavg,yavg,zavg={},{},{},"""
        """Range (x,y,z)={},{},{}""".format(self.quantity, self.modifier_func, self.get_minmaxtime(),
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
                               self.quantity, self.modifier_func, self.get_minmaxtime()))
        else:
            with open(filename, "wb") as file:
                file.write(bytearray("""# Slice diagnostic for {} (with {} applied), t={},\n"""
                                     """# xavg,yavg,zavg={},{},{}, Range (x,y,z)={},{},
                                     {}\n""".format(
                        self.quantity, self.modifier_func, self.get_minmaxtime(),
                        self.diagspace.xavg, self.diagspace.yavg, self.diagspace.zavg,
                        self.diagspace.xslice, self.diagspace.yslice, self.diagspace.zslice),
                        encoding="utf-8"))
                for tind, time in enumerate(self.timearray):
                    np.savetxt(file, np.atleast_1d(np.squeeze(self.dataarray[tind])), fmt='%.9e',
                               header="{}".format(self.timearray[tind]))
                    file.write(b"\n")
