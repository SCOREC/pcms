#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .diagnostic import Diagnostic
from .baseplot import Plotting
from .diagspace import DiagSpace
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import warnings
import utils.averages as averages
from tkinter import END


class diag_contours(Diagnostic):
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.name = 'Contours'
        self.tab = 'Basic'

        self.help_txt = """Plot contours on xy and xz planes
                        \n
                        \nQuantity : quantity to plot 
                        \nx ind : index of radial position (def. all)
                        \ny ind : index of binormal position (def. 0)
                        \nz ind : index of parallel position (def. outboard midplane)
                        \nspec: which species to plot (def. all if spec dependent)
                        \nFourier x : spectral in x direction (def. False)
                        \nFourier y : spectral in y direction (def. False)
                        \nSave h5 : save hdf5 file (def. False)
                        \nt avg : time averaged (def. False)"""

        self.avail_vars = avail_vars
        self.specnames = specnames

        self.opts = {"quant": {'tag': "Quant", 'values': self.ListAvailVars(['mom', 'field'])},
                     "spec": {'tag': "Spec", 'values': self.ListSpecnames()},
                     "xind": {'tag': "x ind", 'values': None},
                     "yind": {'tag': "y ind", 'values': None},
                     "zind": {'tag': "z ind", 'values': None},
                     "f_x": {'tag': "Fourier x", 'values': [False, True]},
                     "f_y": {'tag': "Fourier y", 'values': [False, True]},
                     'save_h5': {'tag': "Save h5", 'values': [True, False]},
                     't_avg': {'tag': "t avg", 'values': [False, True]}, }

        self.set_defaults()

        self.options_gui = Diagnostic.Options_GUI()

    def set_options(self, run_data, species):
        if self.opts['quant']['value'] is None:
            raise RuntimeError("No quantities given for contours")

        """ parallel position, unset mean outboard midplane, -1 all """
        zind = self.opts['zind']['value']
        if not zind:
            zind = int(run_data.pnt.nz0/2)
            zrange = (zind, zind + 1)
        elif zind == -1:
            zind = None
            zrange = (None, None)
        else:
            zrange = (zind, zind + 1)

        """ radial position, unset or -1 means all """
        xind = self.opts['xind']['value']
        if not xind or xind == -1:
            xind = None
            xrange = (None, None)
        else:
            xrange = (xind, xind + 1)

        """ binormal position, unset means zero, -1 means all  """
        yind = self.opts['yind']['value']

        if not yind or yind == -1:
            yind = None
            yrange = (None, None)
        else:
            yrange = (yind, yind + 1)

        """ Fourier in x """
        x_fourier = self.opts['f_x']['value']

        """ Fourier in y """
        y_fourier = self.opts['f_y']['value']

        """ diagspace is unique for each call, so it can be here
            it should however be outside, since we can construct all grids
            and have them in common for all diagnostics and select here only
            the indexes we weant to plot"""
        self.diagspace = DiagSpace(run_data.spatialgrid, x_fourier, y_fourier, False, xrange,
                                   yrange, zrange, False, False, False)

        self.geom = run_data.geometry
        self.run_data = run_data

        """ This is special since we dont save stuff by default, so just do what we need"""
        self.Get_Spec_From_Opts()

        """ toggle the reader """
        self.Get_NeededVars()

        """ cast the output as a dictionary of lists to ease the rest"""
        #        TODO there is probably a best structure for this....
        self.contours = {}
        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():
                if self.needed_vars[file][quant]:
                    if file == 'field':
                        self.contours[quant] = []
                    else:
                        for spec in self.specnames:
                            self.contours[quant + '#' + spec] = []

        return self.needed_vars

    def execute(self, data, step):
        """ not much to do other than appending data """

        """ loop over files """
        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():
                """ loop over quaitites in that file"""
                if self.needed_vars[file][quant]:
                    if file == 'field':
                        """ no species dependency """
                        data_in = getattr(getattr(data, file), quant)(step.time,
                                                                      getattr(step, file))
                        data_in = data._apply_fouriertransforms(self.diagspace, data_in)
                        self.contours[quant].append(data_in[self.diagspace.diagslice])
                    else:
                        """ spec dependent"""
                        for spec in self.specnames:
                            data_in = getattr(getattr(data, file + '_' + spec), quant)(step.time,
                                                                                       getattr(step,
                                                                                               file))

                            # dont know if we should use this
                            data_in = averages.av3d_by_switch(self.diagspace.xavg,
                                                              self.diagspace.yavg,
                                                              self.diagspace.zavg)(data_in,
                                                                                   self.geom)
                            data_in = data._apply_fouriertransforms(self.diagspace, data_in)
                            self.contours[quant + '#' + spec].append(
                                    data_in[self.diagspace.diagslice])

    def plot(self, time_requested, output=None):
        """ Not much to explain """

        def round_to_n(x, n):
            # Helper function to round to n significant digits
            return np.round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))

        def plot_a_quant(ax, x, y, f, is_f, x_lbl, y_lbl, ttl):

            absval = int(round_to_n(np.max(np.abs(f) if is_f else f), 2))
            cm1 = ax.pcolor(x, y, np.abs(np.squeeze(f)).T if is_f else np.squeeze(f).T)
            #            cm1 = ax.contourf(x, y, np.abs(np.squeeze(f)).T if is_f else np.squeeze(
            #            f).T,
            #                            100, cmap=self.plotbase.cmap_bidirect)
            ax.set_rasterization_zorder(z=-10)
            ax.tick_params(axis='both', which='major', labelsize=self.plotbase.xyfs)
            ax.set_xlabel(x_lbl)
            ax.set_ylabel(y_lbl)
            ax.set_title(ttl)
            fig.colorbar(cm1)

        self.plotbase = Plotting()
        self.plotbase.titles.update(
                {"phi": r"$\phi$", "A_par": r"$A_{//}$", "B_par": r"$B_{//}$", "dens": r"$n$",
                 "T_par": r"T_{//}", "T_perp": r"T_{\perp}", "u_par": r"u_{//}", "q_par": r"q_{//}",
                 "q_perp": r"q_{\perp}"})

        if self.opts['f_x']['value']:
            x_lbl = r"$k_x\rho_{ref}$"
            x = self.run_data.spatialgrid.kx_fftorder
        else:
            x_lbl = r"$x/\rho_{ref}$"
            x = self.run_data.spatialgrid.x

        if self.opts['f_y']['value']:
            y_lbl = r"$k_y\rho_{ref}$"
            y = self.run_data.spatialgrid.ky
        else:
            y_lbl = r"$y/\rho_{ref}$"
            y = self.run_data.spatialgrid.y

        is_f = self.opts['f_x']['value'] or self.opts['f_y']['value']

        if self.opts['t_avg']['value']:
            """ loop over file"""
            for quant in self.contours.keys():
                fig = plt.figure()
                ind_str = quant.find('#')
                ttl = self.plotbase.titles[quant] if ind_str == -1 else self.plotbase.titles[quant[
                                                                                             0:ind_str]] + " " + quant[
                                                                                                                 ind_str + 1:]
                ax = fig.add_subplot(111)
                plot_a_quant(ax, x, y, averages.mytrapz(self.contours[quant], time_requested), is_f,
                             x_lbl, y_lbl, ttl)
                fig.show()

        elif len(time_requested) < 12:
            """ in this case we plot on a single figure for each quantity"""
            nr = int(np.ceil(np.sqrt(len(time_requested))))
            nc = int(np.ceil(len(time_requested)/nr))
            """ loop over file"""
            for quant in self.contours.keys():
                fig = plt.figure()

                ind = quant.find('#')
                if ind == -1:
                    ttl = self.plotbase.titles[quant]
                else:
                    fig.clear()
                    fig = plt.figure()
                    ttl = self.plotbase.titles[quant[0:ind]] + " " + quant[ind + 1:]
                for i_t, time in enumerate(time_requested):
                    ax = fig.add_subplot(nr, nc, i_t + 1)
                    plot_a_quant(ax, x, y, self.contours[quant][i_t], is_f, x_lbl, y_lbl,
                                 ttl + " @ t={:.3f}".format(time))
                # with tight_layour the figures dont look good.
                # fig.tight_layout()
                fig.show()

        else:
            warnings.warn("Too many timesteps selected", RuntimeWarning)
            if output:
                output.info_txt.insert(END, "Too many timesteps selected\n")
                output.info_txt.see(END)
