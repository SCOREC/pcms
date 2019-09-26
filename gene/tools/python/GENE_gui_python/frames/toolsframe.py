#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import LabelFrame, Button, BOTTOM, BOTH, Toplevel, TOP, END, Scrollbar, Text, RIGHT, \
    LEFT, Y
import tkinter as tk
import matplotlib
from data.nrgdata import gluenrgdata
import numpy as np

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import math

LARGE_FONT = ("Helvetica", 12)

"""   Class  containing the basic tools of a given simualtion

plots geomety, profiles, time traces, and so on"""


class ToolsFrame:

    def __init__(self, parent, gui, sim):

        self.plotframe = LabelFrame(parent, text="Tools")
        self.plotframe.grid(row=0, column=0, sticky="nswe")

        self.tracesbtn = Button(self.plotframe, text="Plot Traces",
                                command=lambda: self.plot_traces(parent, gui, sim))
        self.tracesbtn.grid(row=1, column=0, sticky="nswe", padx=10, pady=10)

        self.geometrybtn = Button(self.plotframe, text="Plot Geometry",
                                  command=lambda: self.plot_geometry(parent, gui, sim))
        self.geometrybtn.grid(row=1, column=1, sticky="nswe", padx=10, pady=10)

        self.profs_btn = Button(self.plotframe, text="Plot Profiles",
                                command=lambda: self.plot_profiles(parent, gui, sim))
        self.profs_btn.grid(row=1, column=2, sticky="nswe", padx=10, pady=10)

        self.infobtn = Button(self.plotframe, text="Informations",
                              command=lambda: self.print_info(parent, gui, sim))
        self.infobtn.grid(row=1, column=3, sticky="nswe", padx=10, pady=10)

        self.refvalues = Button(self.plotframe, text="Reference Values",
                                command=lambda: self.print_refvalues(parent, gui, sim))
        self.refvalues.grid(row=2, column=0, sticky="nswe", padx=10, pady=10)

    def plot_traces(self, parent, gui, sim):
        """ create one figure only where we plot the fluxes if is nonlinear
            or fluctuation traces if is linear. For more stuff there is the nrg tab """
        if not sim.extensions:
            gui.status.info_txt.insert(END, "... maybe you should select a run first?\n")
        else:
            window = Toplevel(parent)
            tk.Tk.wm_title(window, "Time traces")

            # plot based on first ezxtension only
            # TODO  what about a scan?
            if sim.run[0].nonlinear:
                f = self.plot_fluxes(sim)
            else:
                f = self.plot_fluctuations(sim)

            show_this(window, f)

    def plot_fluxes(self, sim):

        n_rows = 3
        n_cols = sim.run[0].pnt.n_spec

        data = gluenrgdata(sim.nrgdata)

        figure = Figure(figsize=(5, 5), dpi=100)
        for i in range(0, n_cols):
            # heat flux
            a = figure.add_subplot(n_rows, n_cols, i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 6])
            leg = [r"$Q_{es}$"]
            if sim.run[0].electromagnetic:
                a.plot(data.timearray, data.dataarray[:, i, 7])
                leg.append(r"$Q_{em}$")
            a.legend(leg)
            a.set_title(sim.run[0].specnames[i], fontsize=20)

            # particle flux
            a = figure.add_subplot(n_rows, n_cols, n_cols + i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 4])
            leg = [r"$\Gamma_{es}$"]
            if sim.run[0].electromagnetic:
                a.plot(data.timearray, data.dataarray[:, i, 5])
                leg.append(r"$\Gamma_{em}$")
            a.legend(leg)

            # momentum flux
            a = figure.add_subplot(n_rows, n_cols, n_cols*2 + i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 8])
            leg = [r"$\Pi_{es}$"]
            if sim.run[0].electromagnetic:
                a.plot(data.timearray, data.dataarray[:, i, 9])
                leg.append(r"$\Pi_{em}$")
            a.legend(leg)
            a.set_xlabel(r'$t\: [L_\mathrm{ref}/c_\mathrm{ref}]$', fontsize=15)

        return figure

    def plot_fluctuations(self, sim):

        n_rows = 4
        n_cols = sim.run[0].pnt.n_spec

        data = gluenrgdata(sim.nrgdata)

        figure = Figure(figsize=(5, 5), dpi=100)
        for i in range(0, n_cols):
            a = figure.add_subplot(n_rows, n_cols, i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 0])
            a.legend([r"$|n^{2}|$"])
            a.set_title(sim.run[0].specnames[i], fontsize=20)

            a = figure.add_subplot(n_rows, n_cols, n_cols + i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 1])
            a.legend([r"$|u_{//}^{2}|$"])

            a = figure.add_subplot(n_rows, n_cols, n_cols*2 + i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 2])
            a.legend([r"$|T_{\parallel}|^{2}$"])

            a = figure.add_subplot(n_rows, n_cols, n_cols*3 + i + 1)
            a.plot(data.timearray, data.dataarray[:, i, 3])
            a.legend([r"$|T_{\perp}|^{2}$"])
            a.set_xlabel(r'$t\: [L_\mathrm{ref}/c_\mathrm{ref}]$', fontsize=15)

        return figure

    def plot_geometry(self, parent, gui, sim):
        """ flux tube: one plot with metric as function of z
            x global: plot as a function of x and z
            xy global for the moment as function of z only"""
        if not sim.extensions:
            gui.status.info_txt.insert(END, "... maybe you should select a run first?\n")
        else:
            window = Toplevel(parent)
            tk.Tk.wm_title(window, "Magnetic geometry")

            # plot based on first extension only  # TODO  what about a scan?
        if sim.run[0].x_local:
            if sim.run[0].y_local:
                f = self.plot_geom_local(sim)
            else:
                print("need to be implemented for yglobal")

        else:
            if sim.run[0].y_local:
                f = self.plot_geom_global(sim)
            else:
                f = self.plot_geom_xyglobal(sim)

        show_this(window, f)

    def plot_geom_local(self, sim):

        n_rows = 2
        n_cols = 6

        self.names = {'gxx': r"$g^{xx}$", 'gxy': r"$g^{xy}$", 'gxz': r"$g^{xz}$",
                      'gyy': r"$g^{yy}$", 'gyz': r"$g^{yz}$", 'gzz': r"$g^{zz}$",
                      'dBdx': r"$dB/dx$", 'dBdy': r"$dB/dy$", 'dBdz': r"$dB/dz$",
                      'jacobian': r"$Jacobian$", 'dxdR': r"$dx/dR$", 'dxdZ': r"$dx/dZ$", }

        figure = Figure(figsize=(5, 5), dpi=100)
        i_plt = 0

        for fld, lbl in self.names.items():
            i_plt += 1
            a = figure.add_subplot(n_rows, n_cols, i_plt)
            a.plot(sim.run[0].spatialgrid.z/np.pi, getattr(sim.run[0].geometry, fld))
            a.set_title(lbl, fontsize=8)
            a.set_xlabel(r'$z/\pi$', fontsize=8)

        return figure

    def plot_geom_global(self, sim):

        n_rows = 2
        n_cols = 6

        self.names = {'gxx': r"$g^{xx}$", 'gxy': r"$g^{xy}$", 'gxz': r"$g^{xz}$",
                      'gyy': r"$g^{yy}$", 'gyz': r"$g^{yz}$", 'gzz': r"$g^{zz}$",
                      'dBdx': r"$dB/dx$", 'dBdy': r"$dB/dy$", 'dBdz': r"$dB/dz$",
                      'jacobian': r"$Jacobian$", 'dxdR': r"$dx/dR$", 'dxdZ': r"$dx/dZ$", }

        figure = Figure(figsize=(5, 5), dpi=100)
        i_plt = 0

        for fld, lbl in self.names.items():
            i_plt += 1
            a = figure.add_subplot(n_rows, n_cols, i_plt, projection='3d')
            X, Z = np.meshgrid(sim.grids[0].x, sim.grids[0].z/np.pi)
            a.plot_surface(Z, X, getattr(sim.geomdata[0], fld)[:, :])
            a.set_title(lbl, fontsize=8)
            a.set_xlabel(r'$z/\pi$', fontsize=8)

        return figure

    def plot_geom_xyglobal(self, sim):

        n_rows = 2
        n_cols = 5

        self.names = {'gxx': r"$g^{xx}$", 'gxy': r"$g^{xy}$", 'gxz': r"$g^{xz}$",
                      'gyy': r"$g^{yy}$", 'gyz': r"$g^{yz}$", 'gzz': r"$g^{zz}$",
                      'dBdx': r"$dB/dx$", 'dBdy': r"$dB/dy$", 'dBdz': r"$dB/dz$",
                      'jacobian': r"$Jacobian$", }

        figure = Figure(figsize=(5, 5), dpi=100)
        i_plt = 0

        for fld, lbl in self.names.items():
            i_plt += 1
            a = figure.add_subplot(n_rows, n_cols, i_plt)
            a.plot(sim.grids[0].z/np.pi, getattr(sim.geomdata[0], fld)[:, 0, 0])
            a.set_title(lbl, fontsize=8)
            a.set_xlabel(r'$z/\pi$', fontsize=8)

        return figure

    def plot_profiles(self, parent, gui, sim):

        if not sim.run[0].x_local:
            n_rows = 4
            n_cols = sim.run[0].pnt.n_spec

            if not sim.extensions:
                gui.status.info_txt.insert(END, "... maybe you should select a run first?\n")
            else:
                window = Toplevel(parent)
                tk.Tk.wm_title(window, "Profiles")

            figure = Figure(figsize=(5, 5), dpi=100)

            for i in range(0, n_cols):
                a = figure.add_subplot(n_rows, n_cols, i + 1)
                a.plot(sim.run[0].profdata.xs, sim.run[0].profdata.T0s[:, i])
                a.legend([r"$T(x)$"])
                a.set_title(sim.run[0].specnames[i], fontsize=20)

                a = figure.add_subplot(n_rows, n_cols, n_cols + i + 1)
                a.plot(sim.run[0].profdata.xs, sim.run[0].profdata.omt0s[:, i])
                a.legend([r"$T'(x)$"])
                a.set_title(sim.run[0].specnames[i], fontsize=20)

                a = figure.add_subplot(n_rows, n_cols, n_cols*2 + i + 1)
                a.plot(sim.run[0].profdata.xs, sim.run[0].profdata.n0s[:, i])
                a.legend([r"$n(x)$"])
                a.set_title(sim.run[0].specnames[i], fontsize=20)

                a = figure.add_subplot(n_rows, n_cols, n_cols*3 + i + 1)
                a.plot(sim.run[0].profdata.xs, sim.run[0].profdata.omn0s[:, i])
                a.legend([r"$n'(x)$"])
                a.set_title(sim.run[0].specnames[i], fontsize=20)

            show_this(window, figure)
        else:
            gui.status.info_txt.insert(END, "You can only plot profiles for global simulations\n")

    def print_info(self, parent, gui, sim):

        if not sim.extensions:
            gui.status.info_txt.insert(END, "... maybe you should select a run first?\n")
        else:
            window = Toplevel(parent)
            tk.Tk.wm_title(window, "Parameters")

        S = Scrollbar(window)
        T = Text(window, font=LARGE_FONT)
        S.pack(side=RIGHT, fill=Y)
        T.pack(side=LEFT, fill=Y)
        S.config(command=T.yview)
        T.config(yscrollcommand=S.set)

        for key, value in sim.run[0].pars.items():

            T.insert(END, ": ".join((str(key), str(value))) + '\n')

    def print_refvalues(self, parent, gui, sim):

        if not sim.extensions:
            gui.status.info_txt.insert(END, "... maybe you should select a run first?\n")
        else:
            window = Toplevel(parent)
            tk.Tk.wm_title(window, "Reference Values")

        S = Scrollbar(window)
        T = Text(window, font=LARGE_FONT)
        S.pack(side=RIGHT, fill=Y)
        T.pack(side=LEFT, fill=Y)
        S.config(command=T.yview)
        T.config(yscrollcommand=S.set)

        mp = 1.672621637E-27
        e = 1.6021766208E-19

        bref = float(sim.run[0].pars["Bref"])*1e4  # (in Gauss)
        nref = float(sim.run[0].pars["nref"])*1e13
        tref = float(sim.run[0].pars["Tref"])*1e3  # (in ev)
        lref = float(sim.run[0].pars["Lref"])  # (in m)
        mref = float(sim.run[0].pars["mref"])  # in proton masses
        qref = 1.0

        lnlambda = 24. - math.log(math.sqrt(nref)/tref)
        beta = 4.03e-11*nref*tref/bref ** 2  # electron beta
        coll = 2.3031e-14*nref/tref ** 2*lref*lnlambda*1e2  # assume Tref=Te

        qref = 1.0
        cref = 9.79e3*math.sqrt(tref/mref)
        gyrofreq = 9.58e3*bref/mref*qref
        rhoref = cref/gyrofreq
        inv_rhostar = lref/rhoref

        # Reference value:
        T.insert(END, "Reference values:\n")
        T.insert(END, "\n")

        T.insert(END, ": ".join(("Bref", str(sim.run[0].pars["Bref"]))) + ' T \n')
        T.insert(END, ": ".join(("nref", str(sim.run[0].pars["nref"]))) + ' e+19 m^3 \n')
        T.insert(END, ": ".join(("Tref", str(tref))) + ' eV \n')
        T.insert(END, ": ".join(("Lref", str(lref))) + ' m \n')
        T.insert(END, ": ".join(("mref", str(mp*mref))) + ' kg \n')
        T.insert(END, ": ".join(("Qref", str(e*qref))) + ' C \n')
        T.insert(END, ": ".join(("Cref", "%g"%cref)) + ' m/s \n')
        T.insert(END, ": ".join(("rhoref", "%g"%rhoref)) + ' m \n')
        T.insert(END, ": ".join(("1/rhostar", "%g"%inv_rhostar)) + ' \n')
        T.insert(END, ": ".join(("beta_ref", "%g"%beta)) + ' \n')
        T.insert(END, ": ".join(("Coll", "%g"%coll)) + ' \n')

        # Derive values:

        for spec in sim.run[0].specnames:

            if (sim.run[0].pars['charge' + spec] < 0):  # for electrons
                ne = (sim.run[0].pars['dens' + spec])*1e13
                te = (sim.run[0].pars['temp' + spec])*1e3
            else:
                if (sim.run[0].pars['charge' + spec] == 1):  # for main ion species
                    ti = (sim.run[0].pars['temp' + spec])*1e3
                    qi = sim.run[0].pars['charge' + spec]
                    mi = sim.run[0].pars['mass' + spec]

        try:
            tau = ti/te
        except:
            # there is only one species, so we put default values
            qi = 1
            ne = 1e13
            te = 1e3
            tau = 1
            mi = 1

        c = 2.99792458e10

        T.insert(END, "\n")
        T.insert(END, "Ions:\n")
        T.insert(END, "\n")

        thermalvel = 9.79e3*math.sqrt(tau*te/mi)
        gyrofreq = 9.58e3*bref/mi*qi
        gyrorad = thermalvel/gyrofreq
        plasmafreq = 1.32e3*qi*math.sqrt(ne/mi)
        skindepth = c/plasmafreq/1e2  # 2.28e7/q*math.sqrt(mi/ne)
        collisionality = 4.8e-8*qi ** 4/math.sqrt(mi)*ne*lnlambda/(tau*te) ** (1.5)

        T.insert(END, ": ".join(("thermal velocity", "%g"%thermalvel)) + ' m/s \n')
        T.insert(END, ": ".join(("gyrofrequency", "%g"%gyrofreq)) + ' rad/s \n')
        T.insert(END, ": ".join(("gyroradious", "%g"%gyrorad)) + ' m \n')
        T.insert(END, ": ".join(("plasma frequency", "%g"%plasmafreq)) + ' rad/s \n')
        T.insert(END, ": ".join(("Skin depth", "%g"%skindepth)) + ' m \n')
        T.insert(END, ": ".join(("Collisionality", "%g"%collisionality)) + ' 1/s \n')

        T.insert(END, "\n")
        T.insert(END, "Electrons:\n")
        T.insert(END, "\n")

        thermalvel_e = 4.19e5*math.sqrt(te)
        gyrofreq_e = 1.76e7*bref
        gyrorad_e = thermalvel_e/gyrofreq_e
        plasmafreq_e = 5.64e4*math.sqrt(ne)
        skindepth_e = c/plasmafreq_e/1e2  # 5.31e5/math.sqrt(ne)
        collisionality_e = 2.91e-6*ne*lnlambda/te ** (1.5)

        T.insert(END, ": ".join(("thermal velocity", "%g"%thermalvel_e)) + ' m/s \n')
        T.insert(END, ": ".join(("gyrofrequency", "%g"%gyrofreq_e)) + ' rad/s \n')
        T.insert(END, ": ".join(("gyroradious", "%g"%gyrorad_e)) + ' m \n')
        T.insert(END, ": ".join(("plasma frequency", "%g"%plasmafreq_e)) + ' rad/s \n')
        T.insert(END, ": ".join(("Skin depth", "%g"%skindepth_e)) + ' m \n')
        T.insert(END, ": ".join(("Collisionality", "%g"%collisionality_e)) + ' 1/s \n')

        T.insert(END, "\n")
        T.insert(END, "Parameters:\n")
        T.insert(END, "\n")

        debye = 7.43e2*math.sqrt(tref/nref)/1e2  # uses electron temperature
        alfvenvel = 2.18e11*bref/math.sqrt(mi*nref)/1e2
        n_Debye = 1.72e9*tref ** (1.5)/math.sqrt(nref)  # uses electron temperature
        beta_i = 4.03e-11*nref*tref/bref ** 2*tau  # ion beta
        debye2_norm = (debye/gyrorad) ** 2  # assume Te=Te
        alfvenvel_c = alfvenvel*1e2/c

        T.insert(END, ": ".join(("Debye lenth", "%g"%debye)) + ' m \n')
        T.insert(END, ": ".join(("n_Debye", "%g"%n_Debye)) + ' \n')
        T.insert(END, ": ".join(("Alfven velocity", "%g"%alfvenvel)) + ' m/s \n')
        T.insert(END, ": ".join(("beta_i", "%g"%beta_i)) + ' \n')
        T.insert(END, ": ".join(("debye2", "%g"%debye2_norm)) + ' \n')
        T.insert(END, ": ".join(("v_A/c", "%g"%alfvenvel_c)) + ' \n')

        T.insert(END, "\n")
        T.insert(END, "GK ordering for ions:\n")
        T.insert(END, "\n")

        T.insert(END, ": ".join(("rho_i/Length", "%g"%(gyrorad/lref))) + ' m \n')
        T.insert(END, ": ".join(
                ("Gyrofrequency/Plasma frequency", "%g"%(gyrofreq/plasmafreq))) + ' m \n')
        T.insert(END, ": ".join(
                ("Collisionalisty/Gyrofrequency", "%g"%(collisionality/gyrofreq))) + ' m \n')
        T.insert(END, ": ".join(("lambda_D/rho_i", "%g"%(debye/gyrorad))) + ' m \n')
        T.insert(END, ": ".join(("v_i/c", "%g"%(thermalvel/c*1e2))) + ' m \n')

        T.insert(END, "\n")
        T.insert(END, "GK ordering for electrons:\n")
        T.insert(END, "\n")

        T.insert(END, ": ".join(("rho_e/Length", "%g"%(gyrorad_e/lref))) + ' m \n')
        T.insert(END, ": ".join(
                ("Gyrofrequency/Plasma frequency", "%g"%(gyrofreq_e/plasmafreq_e))) + ' m \n')
        T.insert(END, ": ".join(
                ("Collisionalisty/Gyrofrequency", "%g"%(collisionality_e/gyrofreq_e))) + ' m \n')
        T.insert(END, ": ".join(("lambda_D/rho_i", "%g"%(debye/gyrorad_e))) + ' m \n')
        T.insert(END, ": ".join(("v_i/c", "%g"%(thermalvel_e/c*1e2))) + ' m \n')


class show_this:
    def __init__(self, window, f):
        canvas = FigureCanvasTkAgg(f, window)
        canvas.draw()
        canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

        try:
            from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
            toolbar = NavigationToolbar2Tk(canvas, window)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
        except:
            pass
