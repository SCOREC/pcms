#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .diagnostic import Diagnostic
from .baseplot import Plotting
import matplotlib.pyplot as plt
import numpy as np
import warnings
import utils.averages as averages
from tkinter import END


class diag_amplitude_spectra(Diagnostic):

    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.name = 'Amplitude spectra'
        self.tab = 'Local'

        self.help_txt = """Plot amplitude spectra in kx/ky space,
                        \naveraged over the remaining coordinates and times.
                        \n
                        \nQuantity : quantity to plot 
                        \nspec: which species to plot (def. all if spec dependent)
                        \nSave h5 : save hdf5 file (def. False)
                        """

        self.avail_vars = avail_vars
        self.specnames = specnames

        self.opts = {"quant": {'tag': "Quant", 'values': self.ListAvailVars(['mom', 'field'])},
                     "spec": {'tag': "Spec", 'values': self.ListSpecnames()},
                     'save_h5': {'tag': "Save h5", 'values': [True, False]}, }

        self.set_defaults()

        self.options_gui = Diagnostic.Options_GUI()

    def set_options(self, run_data, specnames):

        if self.opts['quant']['value'] is None:
            raise RuntimeError("No quantities given for contours")

        self.geom = run_data.geometry
        self.run_data = run_data

        """ This is special since we dont save stuff by default, so just do what we need"""
        self.Get_Spec_From_Opts()

        """ toggle the reader """
        self.Get_NeededVars()
        self.amplitude_kx = {}
        self.amplitude_ky = {}

        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():

                if self.needed_vars[file][quant]:

                    if file == 'field':
                        self.amplitude_kx[quant] = []
                        self.amplitude_ky[quant] = []

                    else:
                        for spec in self.specnames:
                            self.amplitude_kx[quant + '#' + spec] = []
                            self.amplitude_ky[quant + '#' + spec] = []

        return self.needed_vars

    def execute(self, data, step):

        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():
                """ loop over quaitites in that file"""
                if self.needed_vars[file][quant]:

                    if file == 'field':
                        """ no species dependency """
                        data_in = getattr(getattr(data, file), quant)(step.time,
                                                                      getattr(step, file))
                        data_in = data_in*np.conj(data_in)
                        self.amplitude_kx[quant].append(
                            averages.flux_spectra_yz_av(data_in, self.geom))
                        self.amplitude_ky[quant].append(
                            averages.flux_spectra_xz_av(data_in, self.geom))


                    else:
                        """ spec dependent"""
                        for spec in self.specnames:
                            data_in = getattr(getattr(data, file + '_' + spec), quant)(step.time,
                                                                                       getattr(step,
                                                                                               file))
                            data_in = data_in*np.conj(data_in)
                            self.amplitude_kx[quant + '#' + spec].append(
                                averages.flux_spectra_yz_av(data_in, self.geom))
                            self.amplitude_ky[quant + '#' + spec].append(
                                averages.flux_spectra_xz_av(data_in, self.geom))

    def plot(self, time_requested, output=None):
        """ For each selected species we have one figure with six subplots.
            Left is vs. kx, right vs. ky; columnwise we plot, log-log, log-lin, lin-in
            Dashed lines are negative values in log-log plot
            Dashed lines are k multiplied values in log-lin plot"""

        if output:
            output.info_txt.insert(END, "Amplitude specta:\n")

        self.plotbase = Plotting()
        self.plotbase.titles.update(
                {"Ges": r"$\Gamma_{es}$", "Qes": r"$Q_{es}$", "Pes": r"$\Pi_{es}$",
                 "Gem": r"$\Gamma_{em}$", "Qem": r"$Q_{em}$", "Pem": r"$\Pi_{em}$"})

        kx = self.run_data.spatialgrid.kx_pos
        ky = self.run_data.spatialgrid.ky

        for quant in self.amplitude_kx.keys():

            fig = plt.figure(figsize=(6, 8))

            ax_loglog_kx = fig.add_subplot(3, 2, 1)
            ax_loglin_kx = fig.add_subplot(3, 2, 3)
            ax_linlin_kx = fig.add_subplot(3, 2, 5)
            ax_loglog_ky = fig.add_subplot(3, 2, 2)
            ax_loglin_ky = fig.add_subplot(3, 2, 4)
            ax_linlin_ky = fig.add_subplot(3, 2, 6)

            amplitude_kx = averages.mytrapz(self.amplitude_kx[quant], time_requested)
            amplitude_ky = averages.mytrapz(self.amplitude_ky[quant], time_requested)

            """ log-log =plots, dhased liliens for negative values"""
            baselogkx, = ax_loglog_kx.plot(kx, amplitude_kx)
            baselogky, = ax_loglog_ky.plot(ky, amplitude_ky)

            """ lin-log plots, nothing fancy"""
            baseloglinkx, = ax_loglin_kx.plot(kx, amplitude_kx)
            baseloglinky, = ax_loglin_ky.plot(ky, amplitude_ky)

            """ lin-lin plots, nothing fancy"""
            ax_linlin_kx.plot(kx, amplitude_kx)
            ax_linlin_ky.plot(ky, amplitude_ky)

            """ set things"""
            ax_loglog_kx.loglog()
            ax_loglog_ky.loglog()

            ax_loglin_kx.set_xscale("log")
            ax_loglin_ky.set_xscale("log")

            """ lin-lin plots, nothing fancy"""
            ax_linlin_kx.set_xlim(left=0)
            ax_linlin_kx.set_xlabel(r"$k_x \rho_{ref}$")
            ax_linlin_ky.set_xlabel(r"$k_y \rho_{ref}$")

            ind = quant.find('#')
            if ind == -1:
                ax_loglog_ky.set_title("{}".format(quant))
                ax_loglog_kx.set_title("{}".format(quant))
            else:
                ax_loglog_ky.set_title("{}".format(quant[0:ind] + " " + quant[ind + 1:]))
                ax_loglog_kx.set_title("{}".format(quant[0:ind] + " " + quant[ind + 1:]))
            fig.tight_layout()
            fig.show()
