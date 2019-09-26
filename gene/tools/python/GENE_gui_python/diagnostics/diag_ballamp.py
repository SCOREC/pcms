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


class diag_ballamp(Diagnostic):
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.name = 'Ballooning representation'
        self.tab = 'Local'

        self.help_txt = """Plot the ballooning structure of fields
                        \n
                        \nQuant: which quantity to plot (def. fields)
                        \nSpec: which species to plot (def. all if spec dependent)
                        \nPlot real/imag : plot real and imaginary part (def. True)
                        \nNormalize :  normalize with max amplitude (def. True)
                        \nky ind : index of the mode to plot (def. first non zero)
                        \nt avg : time averaged values (def. False)
                        \nSave h5 : save hdf5 file (def. False)"""

        self.avail_vars = avail_vars
        self.specnames = specnames

        self.opts = {'quant': {'tag': "Quant", 'values': self.ListAvailVars(['mom', 'field'])},
                     "spec": {'tag': "Spec", 'values': self.ListSpecnames()},
                     're_im': {'tag': "Plot real/imag", 'values': [True, False]},
                     'norm': {'tag': "Normalize", 'values': [True, False]},
                     't_avg': {'tag': "t avg", 'values': [False, True]},
                     'kyind': {'tag': "ky ind", 'values': None},
                     'save_h5': {'tag': "Save h5", 'values': [True, False]}, }

        self.set_defaults()

        self.options_gui = Diagnostic.Options_GUI()

    def set_options(self, run_data, species):
        if self.opts['quant']['value'] is None:
            raise RuntimeError("No quantities given for contours")

        if not run_data.x_local or not run_data.y_local:
            raise NotImplementedError("Not yet implemented for nonlocal simulations")

        self.run_data = run_data
        self.species = species

        self.Get_Spec_From_Opts()

        if self.opts['kyind']['value'] is None:
            self.kyind = (1 if self.run_data.pnt.nonlinear else 0)  # default: 1st finite ky
        else:
            self.kyind = self.opts['kyind']['value']

        """ we can set the phase factor here once and for all at intialization"""
        Cyq0_x0 = 1  # self.run_data.geometry.Cy * self.run_data.pnt.q0 / self.run_data.pnt.x0
        sign_shear = (-1 if self.run_data.pnt.shat < 0 else 1)
        nexc_sign = sign_shear*self.run_data.pnt.sign_Ip_CW*self.run_data.pnt.sign_Bt_CW

        if self.run_data.pnt.adapt_lx:
            nexc = 1
        else:
            nexc = int(np.round(self.run_data.pnt.lx*self.run_data.pnt.n_pol*np.absolute(
                self.run_data.pnt.shat)*self.run_data.pnt.kymin*np.absolute(Cyq0_x0)))
        nexc *= nexc_sign
        self.nexcj = nexc*(self.run_data.pnt.ky0_ind + self.kyind)
        self.nconn = int(int(int(self.run_data.pnt.nx0 - 1)/2)/np.abs(self.nexcj))*2 + 1
        try:
            self.phasefac = (-1) ** self.nexcj*np.exp(
                2.0*np.pi*1j*self.run_data.pnt.n0_global*self.run_data.pnt.q0*(
                            self.run_data.pnt.ky0_ind + self.kyind))
        except:
            self.phasefac = 1.0

        self.chi_ext = []
        for ncind in range(-int(self.nconn/2), int(self.nconn/2) + 1):
            self.chi_ext = np.concatenate(
                    [self.chi_ext, (self.run_data.spatialgrid.z + ncind*2.0*np.pi)])

        self.diagspace = DiagSpace(self.run_data.spatialgrid, False, False, False, (None, None),
                                   (None, None), (None, None), False, False, False)

        """ toggle the reader """
        self.Get_NeededVars()

        self.ballamps = {}
        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():
                if self.needed_vars[file][quant]:
                    if file == 'field':
                        self.ballamps[quant] = []
                    else:
                        for spec in self.specnames:
                            self.ballamps[quant + '#' + spec] = []

        return self.needed_vars

    def execute(self, data, step):

        def recast_balloon(opts, data_in):
            if opts['norm']['value']:
                normval = np.array(data_in[0, self.kyind, int(self.run_data.pnt.nz0/2)],
                                   dtype=np.complex128)
                normphase = np.angle(normval)
                normval = np.absolute(normval)
            else:
                normphase = 0.0
                normval = 1.0

            balloon = []
            for ncind in range(-int(self.nconn/2), int(self.nconn/2) + 1):
                balloon = np.concatenate([balloon, np.conj(self.phasefac) ** ncind*data_in[
                    ncind*self.nexcj, self.kyind]])
            return balloon*np.exp(-1j*normphase)/normval

        """ loop over files """
        self.only_field = True
        for file in self.needed_vars.keys():
            """ loop over quaitites in that file"""
            for quant in self.needed_vars[file].keys():
                """ loop over quaitites in that file"""
                if self.needed_vars[file][quant]:
                    if file == 'field':
                        """ no species dependency """
                        data_in = getattr(getattr(data, file), quant)(step.time,
                                                                      getattr(step, file))
                        self.ballamps[quant].append(recast_balloon(self.opts, data_in))
                    else:
                        self.only_field = False
                        """ spec dependent"""
                        for spec in self.specnames:
                            data_in = getattr(getattr(data, file + '_' + spec), quant)(step.time,
                                                                                       getattr(step,
                                                                                               file))
                            self.ballamps[quant + '#' + spec].append(
                                recast_balloon(self.opts, data_in))

    def plot(self, time_requested, output=None):
        self.plotbase = Plotting()
        self.plotbase.titles.update(
                {"phi": r"$\phi$", "A_par": r"$A_{//}$", "B_par": r"$B_{//}$", "dens": r"$n$",
                 "T_par": r"T_{//}", "T_perp": r"T_{\perp}", "u_par": r"u_{//}", "q_par": r"q_{//}",
                 "q_perp": r"q_{\perp}"})

        def plot_a_quant(fig, ncols, i_c, x, y, ttl, quant):

            """ lin-lin plots """
            ax_lin = fig.add_subplot(2, ncols, 1 + 2*(i_c - 1))
            ax_lin.plot(x, np.absolute(y), label="|{}|".format(self.plotbase.titles[quant]),
                        color='black')
            if self.opts['re_im']['value']:
                ax_lin.plot(x, np.real(y), label="Re({})".format(self.plotbase.titles[quant]),
                            color='blue', lw=1)
                ax_lin.plot(x, np.imag(y), label="Im({})".format(self.plotbase.titles[quant]),
                            color='red', lw=1)
            ax_lin.set_xlabel(r"$\chi/\pi$", fontsize=self.plotbase.xyfs)
            ax_lin.set_title(ttl, fontsize=self.plotbase.xyfs)
            ax_lin.axhline(y=0, color='black', lw=1, ls='--')
            ax_lin.legend()

            """ log-log plots """
            ax_log = fig.add_subplot(2, ncols, 2 + 2*(i_c - 1))
            ax_log.plot(x, np.absolute(y), label="|{}|".format(self.plotbase.titles[quant]),
                        color='black')
            if self.opts['re_im']['value']:
                ax_log.plot(x, np.real(y), label="Re({})".format(self.plotbase.titles[quant]),
                            color='blue', lw=1)
                ax_log.plot(x, np.imag(y), label="Im({})".format(self.plotbase.titles[quant]),
                            color='red', lw=1)
            ax_log.set_xlabel(r"$\chi/\pi$", fontsize=self.plotbase.xyfs)
            ax_log.axhline(y=0, color='black', lw=1, ls='--')
            ax_log.set_yscale('log')
            ax_log.legend()

        chi_pi = self.chi_ext/np.pi
        kyval = self.run_data.spatialgrid.ky[self.kyind]
        ttl_ky = r"$k_y\rho={:6.3f} $".format(kyval)

        if self.opts['t_avg']['value']:
            """ loop over file"""
            if self.only_field:
                nc = len(self.ballamps.keys())
                fig = plt.figure()
                for i_q, quant in enumerate(self.ballamps.keys()):
                    plot_a_quant(fig, nc, i_q + 1, chi_pi,
                                 averages.mytrapz(self.ballamps[quant], time_requested), ttl_ky,
                                 quant)
                    fig.show()
            else:
                for quant in self.ballamps.keys():
                    fig = plt.figure()
                    #                    ttl = self.plotbase.titles[quant] if quant.find('#')==-1
                    #                    else self.plotbase.titles[quant[0:quant.find(
                    #                    '#')==-1]]+" "+quant[quant.find('#')==-1+1:]
                    if quant.find('#') == -1:
                        ttl = self.plotbase.titles[quant]
                    else:
                        ttl = quant[quant.find('#') + 1:]
                    plot_a_quant(fig, 1, 1, chi_pi,
                                 averages.mytrapz(self.ballamps[quant], time_requested), ttl,
                                 quant[0:quant.find('#')])
                    fig.show()

        elif len(time_requested) < 12:
            """ in this case we plot on a single figure for each quantity"""
            if self.only_field:
                nc = len(self.ballamps.keys())
                for i_t, time in enumerate(time_requested):
                    fig = plt.figure()
                    for i_q, quant in enumerate(self.ballamps.keys()):
                        plot_a_quant(fig, nc, i_q + 1, chi_pi, self.ballamps[quant][i_t].T,
                                     ttl_ky + "@ t={:.3f}".format(time), quant)
                    fig.show()
            else:
                for quant in self.ballamps.keys():
                    ind = quant.find('#')
                    if ind == -1:
                        ttl = self.plotbase.titles[quant]
                        leg = quant
                    else:
                        leg = quant[0:ind]
                        ttl = quant[ind + 1:] + " "
                    for i_t, time in enumerate(time_requested):
                        fig = plt.figure()
                        plot_a_quant(fig, 1, 1, chi_pi, self.ballamps[quant][i_t].T,
                                     ttl + ttl_ky + "@ t={:.3f}".format(time), leg)
                        fig.show()
        #                   fig.tight_layout()

        else:
            warnings.warn("Too many timesteps selected", RuntimeWarning)
            if output:
                output.info_txt.insert(END, "Too many timesteps selected\n")
                output.info_txt.see(END)
