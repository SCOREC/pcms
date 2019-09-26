#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .diagnostic import Diagnostic
from .baseplot import Plotting
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import warnings
import utils.averages as averages
import utils.derivatives as deriv
from tkinter import END


class diag_shearing_rate(Diagnostic):
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.name = 'Shearing rate'
        self.tab = 'Basic'

        self.help_txt = """Analyze zonal flow and shearing rate.
                        \n
                        \nadd MRS: add MRS to radial plots (def. False)
                        \nx0: at whch position one takes traces (def. center)"""

        self.avail_vars = avail_vars
        self.specnames = specnames

        self.opts = {'MRS': {'tag': "add MRS", 'values': [False, True]},
                     'pos': {'tag': "x0", 'values': None},
                     'GAMs': {'tag': "GAM estimate", 'values': [False, True]}}

        self.set_defaults()

        self.options_gui = Diagnostic.Options_GUI()

    def set_options(self, run_data, species):
        if not run_data.y_local:
            raise NotImplementedError("Not yet implemented for y global simulations")

        self.run_data = run_data
        self.species = species
        self.geom = run_data.geometry
        self.dx = self.run_data.spatialgrid.x[1] - self.run_data.spatialgrid.x[0]

        """ using the methods from diagspace and fourier is cumbersome, so by hand"""
        self.nx0 = self.run_data.spatialgrid.nx0

        self.mypos = self.opts['pos']['value']
        if not self.mypos:
            self.my_pos = np.argmin(np.abs(self.run_data.spatialgrid.x))

        self.shearing_rate = []

        """ toggle the reader """
        self.Get_NeededVars(['field'])

        return self.needed_vars

    def execute(self, data, step):

        class ZonalStep:
            def __init__(self, phi_zonal_x, Er_x, vExB_x, omegaExB_x, abs_phi_fs):
                self.phi_zonal_x = phi_zonal_x
                self.Er_x = Er_x
                self.vExB_x = vExB_x
                self.omegaExB_x = omegaExB_x
                self.abs_phi_fs = abs_phi_fs

        data_kxky = data.field.phi(step.time, step.field)

        if self.run_data.x_local:
            """flux-tube version"""
            data_fsavg = np.squeeze(averages.z_av3d(data_kxky[:, 0, :], self.geom))
            phi_zonal_x = self.nx0*np.real_if_close(np.fft.ifft(data_fsavg, axis=-1), tol=1e-14)
            """ Er """
            Er_x = self.nx0*np.real_if_close(
                np.fft.ifft(-1j*self.run_data.spatialgrid.kx_fftorder*data_fsavg, axis=-1),
                tol=1e-14)
            """ ExB velocity """
            vExB_x = -Er_x/self.geom.Cxy;
            """ shearign rate """
            omegaExB_x = -self.nx0*np.real_if_close(
                np.fft.ifft(np.power(self.run_data.spatialgrid.kx_fftorder, 2)*data_fsavg, axis=-1),
                tol=1e-14)/self.geom.Cxy

            self.shearing_rate.append(
                ZonalStep(phi_zonal_x, Er_x, vExB_x, omegaExB_x, np.abs(data_fsavg)))

        elif not self.run_data.x_local and self.run_data.y_local:
            """ global version    
            we compute the radial electric field, then we fsavg it, as in
            Eq. 5.20 of X.Lapillonne thesis. Then divide by C_xy for ExB
            velocity and further derive for shearing rate. q profile is taken
            into account. This is correct for circular plasma only"""

            phi_zonal_x = np.squeeze(averages.z_av3d(data_kxky[:, 0, :], self.geom))

            """ Er """
            Er_x = deriv.compute(self.dx, phi_zonal_x, 1)
            """ ExB velocity """
            vExB_x = -Er_x/self.geom.Cxy;
            """ shearign rate """
            omegaExB_x = self.run_data.spatialgrid.x/self.geometry.q*deriv.compute(
                self.geometry.q*Er_x/self.run_data.spatialgrid.x, 2)/self.dx

            self.shearing_rate.append(
                ZonalStep(phi_zonal_x, Er_x, vExB_x, omegaExB_x, np.abs(data_fsavg)))

    def plot(self, time_requested, output=None):

        def plot_a_map(ax, x, y, f, x_lbl, y_lbl, ttl):
            cm1 = ax.pcolor(x, y, f)
            #            cm1 = ax.contourf(x, y, f,
            #                            100, cmap=self.plotbase.cmap_bidirect)
            ax.set_rasterization_zorder(z=-10)
            ax.set_xlabel(x_lbl, fontsize=self.plotbase.xyfs)
            ax.set_ylabel(y_lbl, fontsize=self.plotbase.xyfs)
            ax.set_title(ttl, fontsize=self.plotbase.xyfs)
            fig.colorbar(cm1)

        self.plotbase = Plotting()

        x_lbl = r'$x/\rho_{ref}$' if self.run_data.x_local else r'x/a'

        if len(time_requested) > 1:
            """ some maps"""
            fig = plt.figure()
            ax = fig.add_subplot(2, 2, 1)
            plot_a_map(ax, time_requested, self.run_data.spatialgrid.x,
                       np.array([x.phi_zonal_x for x in self.shearing_rate]).T,
                       r'$ t c_{ref}/L_{ref}$ ', x_lbl, r'$ \langle\phi\rangle [c_{ref}/L_{ref}]$')

            ax = fig.add_subplot(2, 2, 2)
            plot_a_map(ax, time_requested, self.run_data.spatialgrid.x,
                       np.array([x.Er_x for x in self.shearing_rate]).T, r'$t c_{ref}/L_{ref}$',
                       x_lbl, r'$E_r [eT_{ref}/ (\rho^*_{ref})^2 L_{ref}]$')

            ax = fig.add_subplot(2, 2, 3)
            plot_a_map(ax, time_requested, self.run_data.spatialgrid.x,
                       np.array([x.vExB_x for x in self.shearing_rate]).T, r'$t c_{ref}/L_{ref}$',
                       x_lbl, r'$v_{ExB} [c_{ref} \rho^*_{ref}]$')

            ax = fig.add_subplot(2, 2, 4)
            plot_a_map(ax, time_requested, self.run_data.spatialgrid.x,
                       np.array([x.omegaExB_x for x in self.shearing_rate]).T,
                       r'$t c_{ref}/L_{ref}$', x_lbl, r'$\omega_{ExB} [c_{ref}/L_{ref}]$')
            fig.show()

            """ time traces """
            fig = plt.figure()
            ax = fig.add_subplot(2 + self.run_data.x_local, 1, 1)
            ax.plot(time_requested, np.array([x.vExB_x[self.my_pos] for x in self.shearing_rate]).T)
            ax.set_xlabel(r'$t c_{ref}/L_{ref}$', fontsize=self.plotbase.xyfs)
            ax.set_ylabel(r'$v_{ExB} [c_{ref} \rho^*_{ref}]$', fontsize=self.plotbase.xyfs)

            ax = fig.add_subplot(2 + self.run_data.x_local, 1, 2)
            ax.plot(time_requested,
                    np.array([x.omegaExB_x[self.my_pos] for x in self.shearing_rate]).T)
            ax.set_xlabel(r'$t c_{ref}/L_{ref}$', fontsize=self.plotbase.xyfs)
            ax.set_ylabel(r'$\omega_{ExB} [c_{ref} \rho^*_{ref}]$', fontsize=self.plotbase.xyfs)

            if self.run_data.x_local:
                ax = fig.add_subplot(2 + self.run_data.x_local, 1, 3)
                ax.plot(time_requested, np.array(
                        [np.sqrt(np.mean(np.power(np.abs(x.omegaExB_x), 2))) for x in
                         self.shearing_rate]).T)
                ax.set_xlabel(r'$t c_{ref}/L_{ref}$', fontsize=self.plotbase.xyfs)
                ax.set_ylabel(r'$\sqrt{|\omega_{ExB}|^2} [c_{ref}/L_{ref}]$',
                              fontsize=self.plotbase.xyfs)

                shear_avg = averages.mytrapz(np.array(
                        [np.sqrt(np.mean(np.power(np.abs(x.omegaExB_x), 2))) for x in
                         self.shearing_rate]), time_requested)
                str_out = "ExB shearing rate= {:.3f}".format(shear_avg)
                if output:
                    output.info_txt.insert(END, str_out + "\n")
                    output.info_txt.see(END)

            fig.show()

        """ zonal spectra"""
        if self.run_data.x_local:
            fig = plt.figure()
            ax = fig.add_subplot(3, 1, 1)
            ax.plot(self.run_data.spatialgrid.kx_pos, averages.mytrapz(np.array(
                    [x.abs_phi_fs[0:len(self.run_data.spatialgrid.kx_pos)] for x in
                     self.shearing_rate]), time_requested))
            ax.set_xlabel(r'$k_x \rho_{ref}$', fontsize=self.plotbase.xyfs)
            ax.set_title(r"$|\langle\phi\rangle|$")
            ax.loglog()

            ax = fig.add_subplot(3, 1, 2)
            ax.plot(self.run_data.spatialgrid.kx_pos, averages.mytrapz(np.array(
                    [x.abs_phi_fs[0:len(self.run_data.spatialgrid.kx_pos)] for x in
                     self.shearing_rate]), time_requested))
            ax.set_xlabel(r'$k_x \rho_{ref}$', fontsize=self.plotbase.xyfs)
            ax.set_xscale("log")

            ax = fig.add_subplot(3, 1, 3)
            ax.plot(self.run_data.spatialgrid.kx_pos, averages.mytrapz(np.array(
                    [x.abs_phi_fs[0:len(self.run_data.spatialgrid.kx_pos)] for x in
                     self.shearing_rate]), time_requested))
            ax.set_xlabel(r'$k_x \rho_{ref}$', fontsize=self.plotbase.xyfs)

            fig.show()

        """ radial plots"""
        fig = plt.figure()
        ax = fig.add_subplot(3, 1, 1)
        ax.plot(self.run_data.spatialgrid.x,
                averages.mytrapz(np.array([x.phi_zonal_x for x in self.shearing_rate]),
                                 time_requested))
        ax.set_xlabel(x_lbl, fontsize=self.plotbase.xyfs)
        ax.set_ylabel(r'$v_{ExB} [c_{ref} \rho^*_{ref}]$', fontsize=self.plotbase.xyfs)

        ax = fig.add_subplot(3, 1, 2)
        ax.plot(self.run_data.spatialgrid.x,
                averages.mytrapz(np.array([x.vExB_x for x in self.shearing_rate]), time_requested))
        ax.set_xlabel(x_lbl, fontsize=self.plotbase.xyfs)
        ax.set_ylabel(r'$v_{ExB} [c_{ref} \rho^*_{ref}]$', fontsize=self.plotbase.xyfs)

        ax = fig.add_subplot(3, 1, 3)
        ax.plot(self.run_data.spatialgrid.x,
                averages.mytrapz(np.array([x.omegaExB_x for x in self.shearing_rate]),
                                 time_requested))
        ax.set_xlabel(x_lbl, fontsize=self.plotbase.xyfs)
        ax.set_ylabel(r'$\omega_{ExB} [c_{ref} \rho^*_{ref}]$', fontsize=self.plotbase.xyfs)

        fig.tight_layout()
