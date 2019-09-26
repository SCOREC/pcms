#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from data.base_file import GENE_file
import warnings
import utils.fourier as fourier
import numpy as np


class Data(object):

    def __init__(self):
        self.avail_vars = {}
        pass

    def load_in(self, run_data, extensions=None):
        class AvailableTimes():
            pass

        class TimeStep(object):
            def __init__(self, times, steps, files):
                self.times = times
                self.steps = steps
                self.files = files

        self.avail_times = AvailableTimes()

        self.run_data = run_data

        """ **********************************"""
        """ add other outputs here  once coded"""
        """ **********************************"""

        """ field file """
        if run_data.pnt.istep_field > 0:
            self.field = GENE_file.create('field', run_data, spec=None, extensions=extensions)
            self.avail_vars = {'field': self.field.my_vars()}
            a, b, c = self.field.set_times_and_inds()
            self.avail_times.field = TimeStep(a, b, c)

        """ mom are the same for all species, so we create an object for each species 
            but one entry for times """
        if run_data.pnt.istep_mom > 0:
            for s in run_data.specnames:
                setattr(self, 'mom_' + s,
                        GENE_file.create('mom', run_data, spec=s, extensions=extensions))

            self.avail_vars.update(
                    {'mom': getattr(getattr(self, 'mom_' + run_data.specnames[0]), 'my_vars')()})
            a, b, c = getattr(getattr(self, 'mom_' + run_data.specnames[0]), 'set_times_and_inds')()
            setattr(self.avail_times, 'mom', TimeStep(a, b, c))

    def _apply_fouriertransforms(self, diagspace, var):
        if diagspace.xavg and diagspace.yavg:  # Nothing to do if both are averaged
            return var
        xaxis = -3  # By default the last three dimensions of var are x, y, z.
        yaxis = -2  # This is changed if averages are applied
        zaxis = -1
        if diagspace.yavg:
            xaxis += 1
            yaxis += 1
        if diagspace.zavg:
            xaxis += 1
            yaxis += 1

        if diagspace.x_fourier != self.run_data.x_local and not diagspace.xavg:
            if self.run_data.x_local:
                var = fourier.kx_to_x(var, self.run_data.pnt.nx0, axis=xaxis)
            else:
                var = fourier.x_to_kx(var, self.run_data.pnt.nx0, axis=xaxis)
        if diagspace.y_fourier != self.run_data.y_local and not diagspace.yavg:
            if self.run_data.y_local:
                var = fourier.ky_to_y(var, self.run_data.pnt.nky0, axis=yaxis)
            else:
                warnings.warn("y-global is not well tested", RuntimeWarning)
                var = fourier.y_to_ky(var, self.run_data.pnt.nky0, axis=yaxis)
        if diagspace.z_fourier:
            var = fourier.z_to_kz(var, self.run_data.pnt.nz0, axis=zaxis)
        return var

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
                                     """# xavg,yavg,zavg={},{},{}, Range (x,y,z)={},{},
                                     {}\n""".format(self.varname, self.modifier_func,
                        self.get_minmaxtime(), self.diagspace.xavg, self.diagspace.yavg,
                        self.diagspace.zavg, self.diagspace.xslice, self.diagspace.yslice,
                        self.diagspace.zslice), encoding="utf-8"))
                for tind, time in enumerate(self.timearray):
                    np.savetxt(file, np.atleast_1d(np.squeeze(self.dataarray[tind])), fmt='%.9e',
                               header="{}".format(self.timearray[tind]))
                    file.write(b"\n")
