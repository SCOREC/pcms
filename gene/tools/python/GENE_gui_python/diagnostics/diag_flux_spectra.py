from .diagnostic import Diagnostic
from .baseplot import Plotting
import matplotlib.pyplot as plt
import numpy as np
import warnings
import utils.averages as averages
from tkinter import END


class diag_flux_spectra(Diagnostic):

    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.name = 'Flux spectra'
        self.tab = 'Local'

        self.help_txt = """Plot flux spectra in kx/ky space,
                        \naveraged over the remaining coordinates and times.
                        \nDotted lines in the log/log plot indicate negative values.
                        \n
                        \nspec: which species to plot (def. all if spec dependent)
                        \nSave h5 : save hdf5 file (def. False)
                        """

        self.avail_vars = avail_vars
        self.specnames = specnames

        self.opts = {"spec": {'tag': "Spec", 'values': self.ListSpecnames()},
                     "save_h5": {"tag": "Save h5", "values": [True, False]}}

        self.set_defaults()

        self.options_gui = Diagnostic.Options_GUI()

    def set_options(self, run_data, specnames):

        class Flux_Step():
            def __init__(self, specnames, electromagnetic):
                self.__em__ = electromagnetic
                for spec in specnames:
                    setattr(self, spec, self.__Spec_Spectra__(self.__em__))
                pass

            class __Spec_Spectra__():
                def __init__(self, em):
                    self.Qes = self.__FluxSpectra__()
                    self.Ges = self.__FluxSpectra__()
                    self.Pes = self.__FluxSpectra__()

                    if em:
                        self.Qem = self.__FluxSpectra__()
                        self.Gem = self.__FluxSpectra__()
                        self.Pem = self.__FluxSpectra__()

                class __FluxSpectra__():
                    def __init__(self):
                        self.kx = []
                        self.ky = []

        self.geom = run_data.geometry
        self.run_data = run_data
        self.specnames = self.run_data.specnames

        self.flux_step = Flux_Step(self.specnames, self.run_data.electromagnetic)

        self.Get_NeededVars(['field', 'mom'])

        """ This is special since we dont save stuff by default, so just do what we need"""
        self.Get_Spec_From_Opts()

        return self.needed_vars

    def execute(self, data, step):

        def compute_kxky(flux_spec, var, geometry):
            flux_spec.kx.append(averages.flux_spectra_yz_av(var, geometry))
            flux_spec.ky.append(averages.flux_spectra_xz_av(var, geometry))
            return flux_spec

        vE_x = -1j*self.run_data.spatialgrid.ky[np.newaxis, :, np.newaxis]*data.field.phi(step.time,
                                                                                          step.field)
        if self.run_data.electromagnetic:
            B_x = 1j*self.run_data.spatialgrid.ky[np.newaxis, :, np.newaxis]*data.field.A_par(
                step.time, step.field)

        for spec in self.specnames:

            n0 = self.run_data.pars["dens{}".format(spec)]
            T0 = self.run_data.pars["temp{}".format(spec)]
            mass = self.run_data.pars["mass{}".format(spec)]
            charge = self.run_data.pars["charge{}".format(spec)]

            dens = getattr(getattr(data, 'mom_' + spec), "dens")(step.time, step.mom)
            T_par = getattr(getattr(data, 'mom_' + spec), "T_par")(step.time, step.mom)
            T_perp = getattr(getattr(data, 'mom_' + spec), "T_perp")(step.time, step.mom)
            u_par = getattr(getattr(data, 'mom_' + spec), "u_par")(step.time, step.mom)

            flux_step = getattr(self.flux_step, spec)

            G_es = vE_x*np.conj(dens)*n0/self.geom.Cxy
            flux_step.Ges = compute_kxky(flux_step.Ges, G_es, self.geom)

            Q_es = (vE_x*np.conj(0.5*T_par + T_perp + 3/2*dens)*n0*T0)/self.geom.Cxy
            flux_step.Qes = compute_kxky(flux_step.Qes, Q_es, self.geom)

            P_es = (vE_x*np.conj(u_par)*n0*mass)/self.geom.Cxy
            flux_step.Pes = compute_kxky(flux_step.Pes, P_es, self.geom)

            if self.run_data.electromagnetic:
                Q_par = getattr(getattr(data, 'mom_' + spec), "q_par")(step.time, step.mom)
                Q_perp = getattr(getattr(data, 'mom_' + spec), "q_perp")(step.time, step.mom)

                G_em = B_x*np.conj(u_par)*n0//self.geom.Cxy
                flux_step.Gem = compute_kxky(flux_step.Gem, G_em, self.geom)
                Q_em = B_x*np.conj(Q_par + Q_perp)*n0*T0
                flux_step.Qem = compute_kxky(flux_step.Qem, Q_em, self.geom)

                P_em = B_x*(T_par + dens)*n0*mass//self.geom.Cxy
                flux_step.Pem = compute_kxky(flux_step.Pem, P_em, self.geom)

                if self.run_data.bpar:
                    # todo check if the normalization is correct
                    dBpar_dy_norm = -1j*self.run_data.spatialgrid.ky[np.newaxis, :,
                                        np.newaxis]*data.field.B_par(step.time,
                                                                     step.field)/self.geom.Bfield[
                                                                                 np.newaxis,
                                                                                 np.newaxis,
                                                                                 :]*T0/charge
                    densI1 = getattr(getattr(data, 'mom_' + spec), "densI1")(step.time, step.mom)
                    TparI1 = getattr(getattr(data, 'mom_' + spec), "TparI1")(step.time, step.mom)
                    TppI1 = getattr(getattr(data, 'mom_' + spec), "TppI1")(step.time, step.mom)

                    G_em = G_em + (dBpar_dy_norm*np.conj(densI1))*n0/self.geom.Cxy
                    Q_em = Q_em + (dBpar_dy_norm*np.conj(TparI1 + TppI1))*n0*T0/self.geom.Cxy

    def plot(self, time_requested, output=None):
        """ For each selected species we have one figure with six subplots.
            Left is vs. kx, right vs. ky; columnwise we plot, log-log, log-lin, lin-in
            Dashed lines are negative values in log-log plot
            Dashed lines are k multiplied values in log-lin plot"""

        if output:
            output.info_txt.insert(END, "Flux specta:\n")

        self.plotbase = Plotting()
        self.plotbase.titles.update(
                {"Ges": r"$\Gamma_{es}$", "Qes": r"$Q_{es}$", "Pes": r"$\Pi_{es}$",
                 "Gem": r"$\Gamma_{em}$", "Qem": r"$Q_{em}$", "Pem": r"$\Pi_{em}$"})

        kx = self.run_data.spatialgrid.kx_pos
        ky = self.run_data.spatialgrid.ky

        for spec in self.specnames:
            fig = plt.figure(figsize=(6, 8))

            ax_loglog_kx = fig.add_subplot(3, 2, 1)
            ax_loglin_kx = fig.add_subplot(3, 2, 3)
            ax_linlin_kx = fig.add_subplot(3, 2, 5)
            ax_loglog_ky = fig.add_subplot(3, 2, 2)
            ax_loglin_ky = fig.add_subplot(3, 2, 4)
            ax_linlin_ky = fig.add_subplot(3, 2, 6)

            spec_flux = getattr(self.flux_step, spec)

            for flux in vars(spec_flux).keys():

                flux_kx = averages.mytrapz(getattr(getattr(spec_flux, flux), 'kx'), time_requested)
                flux_ky = averages.mytrapz(getattr(getattr(spec_flux, flux), 'ky'), time_requested)
                # Mask negative flux values for solid lines

                pos_flux_kx = np.ma.masked_where((flux_kx <= 0), flux_kx)
                # Mask zero flux for dashed lines, this takes care of the Nyquist mode in kx
                all_flux_kx = np.ma.masked_where((flux_kx == 0), flux_kx)

                pos_flux_ky = np.ma.masked_where((flux_ky <= 0), flux_ky)
                all_flux_ky = np.ma.masked_where((flux_ky == 0), flux_ky)

                """ log-log =plots, dhased liliens for negative values"""
                baselogkx, = ax_loglog_kx.plot(kx, pos_flux_kx, label=self.plotbase.titles[flux])
                ax_loglog_kx.plot(kx, np.abs(all_flux_kx), ls="--", color=baselogkx.get_color())
                baselogky, = ax_loglog_ky.plot(ky, pos_flux_ky, label=self.plotbase.titles[flux])
                ax_loglog_ky.plot(ky, np.abs(all_flux_ky), ls="--", color=baselogky.get_color())

                """ lin-log plots, nothing fancy"""
                baseloglinkx, = ax_loglin_kx.plot(kx, all_flux_kx, label=self.plotbase.titles[flux])
                ax_loglin_kx.plot(kx, all_flux_kx*kx, ls="--", color=baseloglinkx.get_color())
                baseloglinky, = ax_loglin_ky.plot(ky, all_flux_ky, label=self.plotbase.titles[flux])
                ax_loglin_ky.plot(ky, all_flux_ky*ky, ls="--", color=baseloglinky.get_color())

                """ lin-lin plots, nothing fancy"""
                ax_linlin_kx.plot(kx, all_flux_kx, label=self.plotbase.titles[flux])
                ax_linlin_ky.plot(ky, all_flux_ky, label=self.plotbase.titles[flux])

                str_out = "{} {} = {:.4f} (kx intg.) - {:.4f} (ky intg.)".format(spec, flux,
                                                                                 np.sum(flux_kx),
                                                                                 np.sum(flux_ky))
                if output:
                    output.info_txt.insert(END, str_out + "\n")
                    output.info_txt.see(END)

            """ set things"""
            ax_loglog_kx.loglog()
            ax_loglog_kx.set_xlabel(r"$k_x \rho_{ref}$")

            ax_loglog_ky.loglog()
            ax_loglog_ky.set_xlabel(r"$k_y \rho_{ref}$")

            ax_loglin_kx.set_xscale("log")
            ax_loglin_ky.set_xscale("log")

            """ lin-lin plots, nothing fancy"""
            ax_linlin_kx.set_xlim(left=0)
            ax_linlin_kx.set_xlabel(r"$k_x \rho_{ref}$")

            ax_linlin_ky.set_xlabel(r"$k_y \rho_{ref}$")

            for ax in [ax_loglog_kx, ax_loglin_kx, ax_linlin_kx, ax_loglog_ky, ax_loglin_ky,
                       ax_linlin_ky, ]:
                # ax.set_ylabel(r"$<|A|^2>$")
                ax.legend()
            ax_loglog_ky.set_title("{}".format(spec))
            ax_loglog_kx.set_title("{}".format(spec))
            #            fig.tight_layout()
            fig.show()
