""" Classes to do plots of kx and ky spectra"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
from multiprocessing.dummy import Pool

from pydiag.diagplots.baseplot import Plotting
from pydiag.data.slices import MomFieldSlice
from pydiag.utils.comm import DiagSpace

import pydiag.utils.averages as averages
import pydiag.utils.geom


class AmplitudeSpectra(object):
    """ kx and ky spectra of 3d mom and field data

     Collect the squared absolute field and mom data as kx and ky spectra
     """

    def __init__(self, common, species, rundatafiles, moms=None):
        if not common.y_local:
            raise NotImplementedError("Only implemented for local simulations yet")
        self.cm = common
        self.species = species
        self.diagspace_kx = DiagSpace(common.spatialgrid, x_fourier=True, y_fourier=True, yavg=True,
                                      zavg=True)
        self.diagspace_ky = DiagSpace(common.spatialgrid, x_fourier=True, y_fourier=True, xavg=True,
                                      zavg=True)
        self.momamps_kx = {"phi": None, "dens": None, "tpar": None, "tperp": None, "qpar": None,
                           "qperp": None, "upar": None}
        self.momamps_ky = {"phi": None, "dens": None, "tpar": None, "tperp": None, "qpar": None,
                           "qperp": None, "upar": None}
        if common.electromagnetic:
            self.momamps_kx.update({"apar": None})
            self.momamps_ky.update({"apar": None})
            if common.bpar:
                self.momamps_kx.update({"bpar": None})
                self.momamps_ky.update({"bpar": None})
        if moms:
            self.momamps_kx = {mom: self.momamps_kx[mom] for mom in moms if mom in self.momamps_kx}
            self.momamps_ky = {mom: self.momamps_ky[mom] for mom in moms if mom in self.momamps_ky}

        with Pool(len(self.momamps_kx)) as p:
            momamps_kx_temp = p.map(
                lambda mom: MomFieldSlice(common, mom, species, self.diagspace_kx, rundatafiles),
                self.momamps_kx)
            momamps_ky_temp = p.map(
                lambda mom: MomFieldSlice(common, mom, species, self.diagspace_ky, rundatafiles),
                self.momamps_ky)
        for imom, mom in enumerate(self.momamps_kx):
            self.momamps_kx[mom] = momamps_kx_temp[imom]
            self.momamps_kx[mom].generate_timeseries()
            self.momamps_kx[mom].generate_timeaverage()
            self.momamps_ky[mom] = momamps_ky_temp[imom]
            self.momamps_ky[mom].generate_timeseries()
            self.momamps_ky[mom].generate_timeaverage()

class PerpSpectra(object):
    """ 1d perpendicular spectra of 3d mom and field data

     Collect the squared absolute field and mom data as kperp spectra
     """

    def __init__(self, common, species, rundatafiles):
        if not common.x_local or not common.y_local:
            raise NotImplementedError("Only implemented for local simulations yet")
        self.cm = common
        self.species = species
        self.diagspace_kxky = DiagSpace(common.spatialgrid, x_fourier=True, y_fourier=True, zavg=True)
        self.momamps_kperp = {"phi": None}#, "dens": None, "tpar": None, "tperp": None, "qpar": None,
#                           "qperp": None, "upar": None}
        if common.electromagnetic:
            self.momamps_kperp.update({"apar": None})
            if common.bpar:
                self.momamps_kperp.update({"bpar": None})

        with Pool() as p:
            momamps_kperp_temp = p.map(
                lambda mom: MomFieldSlice(common, mom, species, self.diagspace_kxky, rundatafiles),
                self.momamps_kperp)
        for imom, mom in enumerate(self.momamps_kperp):
            self.momamps_kperp[mom] = momamps_kperp_temp[imom]
            self.momamps_kperp[mom].generate_timeseries()
            self.momamps_kperp[mom].generate_timeaverage()



class FluxSpectra(object):
    """ kx and ky radial flux spectra

    Collect and calculate the kx and ky spectra for the radial fluxes from
    mom files. For further information see diagnostics/doc/spectra.pdf

    """

    def __init__(self, common, species, rundatafiles, with_momentum=True):
        if not common.y_local:
            raise NotImplementedError("Not implemented for y-global simulations yet")
        self.cm = common
        self.species = species
        self.timearray = None
        self.fluxes_kx = {"Gammaes": [], "Qes": []}
        self.fluxes_ky = {"Gammaes": [], "Qes": []}
        if with_momentum:
            self.fluxes_kx.update({"Pies": []})
            self.fluxes_ky.update({"Pies": []})
        if common.electromagnetic:
            self.fluxes_kx.update({"Gammaem": [], "Qem": []})
            self.fluxes_ky.update({"Gammaem": [], "Qem": []})
            if with_momentum:
                self.fluxes_kx.update({"Piem": []})
                self.fluxes_ky.update({"Piem": []})
        self.calc_fluxes(species=species, rundatafiles=rundatafiles, with_momentum=with_momentum)

    def _fetch_moms(self, rundatafiles, species):
        """ Get mom and field slice object for species"""
        diagspace = DiagSpace(self.cm.spatialgrid, x_fourier=True, y_fourier=True)
        momamps = {"phi": None, "dens": None, "tpar": None, "tperp": None, "qpar": None,
                   "qperp": None, "upar": None}
        if self.cm.electromagnetic:
            momamps.update({"apar": None})
            if self.cm.bpar:
                momamps.update({"bpar": None})
        with Pool(len(momamps)) as p:
            momamp_temp = p.map(
                lambda mom: MomFieldSlice(self.cm, mom, species, diagspace, rundatafiles,
                                          modifier_func=lambda dum: dum), momamps)
        for imom, mom in enumerate(momamps):
            momamps[mom] = momamp_temp[imom]
        return momamps

    def calc_fluxes(self, species, rundatafiles, with_momentum):
        """ Calculate the flux spectra"""
        momamps = self._fetch_moms(rundatafiles, species=species)
        ispec = self.cm.specnames.index(species) + 1
        geom = pydiag.utils.geom.Geometry(self.cm)
        momamps["dens"].check_times()
        pos = momamps["dens"].calc_positions()
        self.timearray = np.take(momamps["dens"].timearray, pos)
        for time in self.timearray:
            # vE_x = -c/C_{xy} d phi / dy
            vE_x = - 1j*self.cm.spatialgrid.ky[np.newaxis, :, np.newaxis]*momamps[
                "phi"].generate_slice_attime(time)
            if self.cm.electromagnetic:
                # B_x = d A_par / dy
                B_x = 1j*self.cm.spatialgrid.ky[np.newaxis, :, np.newaxis]*momamps[
                    "apar"].generate_slice_attime(time)
                if self.cm.bpar:
                    # (d B_par / dy) / B0(z)
                    dBpar_dy = 1j*1/self.cm.pnt.Bref*self.cm.spatialgrid.ky[np.newaxis, :,
                                                     np.newaxis]*momamps[
                                   "apar"].generate_slice_attime(time)
                    dBpar_dy /= geom.Bfield[np.newaxis, np.newaxis, :]
                else:
                    dBpar_dy = 0

            def mom_t(momname):
                return momamps[momname].generate_slice_attime(time)

            # Required due to normalisation
            n0 = self.cm.pars["dens{}".format(ispec)]
            mass = self.cm.pars["mass{}".format(ispec)]
            T0 = self.cm.pars["temp{}".format(ispec)]
            # Gamma_es = <n * ve_x>*n0 (n0 due to normalisation)
            # Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
            fluxfactors = {"Gammaes": vE_x*np.conj(mom_t("dens"))*n0, "Qes": vE_x*np.conj(
                    0.5*mom_t("tpar") + mom_t("tperp") + 1.5*mom_t("dens"))*n0*T0}
            if with_momentum:
                warnings.warn("The momentum transport spectra are still incorrect", RuntimeWarning)
                #  P_es = <upar ve_x>
                fluxfactors.update({"Pies": vE_x*mom_t("upar")*n0*mass})
            if self.cm.electromagnetic:
                # Gamma_em = <upar * B_x>
                # Q_em = (qpar + qperp + 5/2 upar) B_x
                # "qpar" and "qperp" from mom file are NOT identical to qpar and qperp
                # but to qpar+1.5upar and qperp+upar
                fluxfactors.update({"Gammaem": B_x*np.conj(mom_t("upar"))*n0,
                                    "Qem": B_x*np.conj(mom_t("qpar") + mom_t("qperp"))*n0*T0})
                if with_momentum:
                    # P_em = (Tpar + n) B_x
                    warnings.warn("The momentum transport spectra are still incorrect",
                                  RuntimeWarning)
                    fluxfactors.update({"Piem": B_x*(mom_t("tpar") + mom_t("dens"))*n0*T0})
                if self.cm.bpar:
                    raise NotImplementedError("bpar part of the fluxes is not implemented yet")
                    fluxfactors["Gammaem"] += 0
                    fluxfactors["Qem"] += 0
            for flux in self.fluxes_kx:
                self.fluxes_kx[flux].append(averages.flux_spectra_yz_av(fluxfactors[flux], geom))
                self.fluxes_ky[flux].append(averages.flux_spectra_xz_av(fluxfactors[flux], geom))
        for flux in self.fluxes_kx:
            self.fluxes_kx[flux] = averages.mytrapz(self.fluxes_kx[flux], self.timearray)
            self.fluxes_ky[flux] = averages.mytrapz(self.fluxes_ky[flux], self.timearray)

    def print_total_fluxes(self):
        """Print the total flux for comparison with nrg diag"""
        print("============================")
        print("Total flux for species {}:".format(self.species))
        for flux in self.fluxes_kx:
            print("{} integrated kx spectrum: {}".format(flux, np.sum(self.fluxes_kx[flux])))
            print("{} integrated ky spectrum: {}".format(flux, np.sum(self.fluxes_ky[flux])))
        print("\n")


class PlotAmplitudeSpectra(Plotting):

    def __init__(self, commonlist, spectraserieslist):
        super().__init__()
        self.spectralist = spectraserieslist
        self.commonlist = commonlist


    def createfigs(self):
        for spectraseries in self.spectralist:
            moms = spectraseries.momamps_kx.keys()
            fig_kx = plt.figure()
            fig_ky = plt.figure()

            ax_log_kx = fig_kx.add_subplot(2, 1, 1)
            ax_lin_kx = fig_kx.add_subplot(2, 1, 2)
            ax_log_ky = fig_ky.add_subplot(2, 1, 1)
            ax_lin_ky = fig_ky.add_subplot(2, 1, 2)
            kx = spectraseries.cm.spatialgrid.kx
            ky = spectraseries.cm.spatialgrid.ky
            for mom in moms:
                ax_log_kx.plot(kx, np.fft.fftshift(spectraseries.momamps_kx[mom].timeaverage),
                               label=self.titles[mom])
                ax_lin_kx.plot(kx, np.fft.fftshift(spectraseries.momamps_kx[mom].timeaverage),
                               label=self.titles[mom])
                ax_log_ky.plot(ky, spectraseries.momamps_ky[mom].timeaverage,
                               label=self.titles[mom])
                ax_lin_ky.plot(ky, spectraseries.momamps_ky[mom].timeaverage,
                               label=self.titles[mom])
            ax_log_kx.loglog()
            ax_log_ky.loglog()
            ax_log_kx.set_xlabel(r"$k_x \rho_{ref}$")
            ax_lin_kx.set_xlabel(r"$k_x \rho_{ref}$")
            ax_log_ky.set_xlabel(r"$k_y \rho_{ref}$")
            ax_lin_ky.set_xlabel(r"$k_y \rho_{ref}$")
            for ax in [ax_log_kx, ax_lin_kx, ax_lin_ky, ax_log_ky]:
                ax.set_ylabel(r"$<|A|^2>$")
                ax.legend()
            ax_log_kx.set_title("{}".format(spectraseries.species))
            ax_log_ky.set_title("{}".format(spectraseries.species))

class PlotPerpSpectra(Plotting):

    def __init__(self, commonlist, spectraserieslist):
        super().__init__()
        self.spectralist = spectraserieslist
        self.commonlist = commonlist
        self.titles.update({"phi": r"$E_\perp$", "apar": r"$B_\perp$", "dens": r"$n$",
                            "tpar": r"$T_\parallel$", "tperp": r"$T_\perp$",
                            "qpar": r"$q_\parallel + 1.5p_0 u_\parallel$",
                            "qperp": r"$q_\perp + p_0 u_\parallel$", "upar": r"$u_\parallel$",
                            "bpar": r"$B_\parallel$", "n00": r"$N_{00}$",
                            "n20": r"$N_{20}$", "n02": r"$N_{02}$"})

    def createfigs(self, moms=None):
        for spectraseries in self.spectralist:
            if not moms:
                moms = spectraseries.momamps_kperp
            fig_kperp = plt.figure()

            ax_log_kperp = fig_kperp.add_subplot(1, 1, 1)

            kx = spectraseries.cm.spatialgrid.kx_fftorder
            ky = spectraseries.cm.spatialgrid.ky
            kperp_sq_2d = kx[:,None]**2+ky[None,:]**2

            #generate kperp array for binning, default number of bins is 1.5*nky0
            nbins=int(np.rint(0.5*len(ky)))
            kperp_bins=np.linspace(0,np.sqrt(np.amax(kperp_sq_2d)),nbins+1)
            #kperp bin centers for plotting (match also array sizes)
            kperp_bc=kperp_bins[:-1]+0.5*np.diff(kperp_bins)[0]

            #logarithmic binning -- disabled by default
            #kperp_bins=np.logspace(np.log10(0.5*ky[1]),np.log10(np.sqrt(np.amax(kperp_sq_2d))),nbins+1)
            #bin "centers" as geometric means of bin boundaries
            #kperp_bc=np.sqrt(kperp_bins[:-1]*kperp_bins[1:])

            for mom in moms:
                if mom in ("phi","apar"):
                    tmp=kperp_sq_2d*spectraseries.momamps_kperp[mom].timeaverage
                else:
                    tmp=spectraseries.momamps_kperp[mom].timeaverage
                #generate 1d perpendicular spectra by binning 2d kperp values
                #with spectral amplitude as the weight
                mom_kperp,kperp_bins=np.histogram(np.sqrt(kperp_sq_2d),bins=kperp_bins,weights=tmp)
                ax_log_kperp.plot(kperp_bc, mom_kperp,label=self.titles[mom])
            xmin,xmax=kperp_bc[0],kperp_bc[np.int(nbins/np.sqrt(2))]
            ax_log_kperp.set_xlim(xmin,xmax)
            kp_mhd=np.linspace(xmin,1,128)
            kp_kin=np.linspace(1,xmax,128)
            ax_log_kperp.plot(kp_mhd, kp_mhd**(-5./3),ls='--',lw=1,label="$k_\perp^{-5/3}$")
            ax_log_kperp.plot(kp_kin, kp_kin**(-2.8),ls='--',lw=1,label="$k_\perp^{-2.8}$")
            #ax_log_kperp.plot(kperp_bc, kperp_bc**(-7./3),ls='--',lw=1,label="$k_\perp^{-7/3}$")
            ax_log_kperp.loglog()
            ax_log_kperp.set_xlabel(r"$k_\perp \rho_{ref}$")
            ax_log_kperp.set_ylabel(r"$\langle|A|^2\rangle$")
            ax_log_kperp.legend()
            ax_log_kperp.set_title("{}".format('Field energy spectra'))
            fig_kperp.tight_layout()


class PlotFluxSpectra(Plotting):

    def __init__(self, commonlist, fluxspecserieslist):
        super().__init__()
        self.fluxspeclist = fluxspecserieslist
        self.commonlist = commonlist
        self.titles.update({"Gammaes": r"$\Gamma_{es}$", "Qes": r"$Q_{es}$", "Pies": r"$\Pi_{es}$",
                            "Gammaem": r"$\Gamma_{em}$", "Qem": r"$Q_{em}$", "Piem": r"$\Pi_{em}$"})

    def createfigs(self):
        for fluxes in self.fluxspeclist:
            fig_kx = plt.figure(figsize=(6, 8))
            fig_ky = plt.figure(figsize=(6, 8))

            ax_log_kx = fig_kx.add_subplot(2, 1, 1)
            ax_lin_kx = fig_kx.add_subplot(2, 1, 2)
            ax_log_ky = fig_ky.add_subplot(2, 1, 1)
            ax_lin_ky = fig_ky.add_subplot(2, 1, 2)
            kx = fluxes.cm.spatialgrid.kx_pos
            ky = fluxes.cm.spatialgrid.ky
            for flux in fluxes.fluxes_kx:
                # Mask negative flux values for solid lines
                pos_flux_kx = np.ma.masked_where((fluxes.fluxes_kx[flux] <= 0),
                                                 fluxes.fluxes_kx[flux])
                # Mask zero flux for dashed lines, this takes care of the Nyquist mode in kx
                all_flux_kx = np.ma.masked_where((fluxes.fluxes_kx[flux] == 0),
                                                 fluxes.fluxes_kx[flux])
                pos_flux_ky = np.ma.masked_where((fluxes.fluxes_ky[flux] <= 0),
                                                 fluxes.fluxes_ky[flux])
                all_flux_ky = np.ma.masked_where((fluxes.fluxes_ky[flux] == 0),
                                                 fluxes.fluxes_ky[flux])
                baselogkx, = ax_log_kx.plot(kx, pos_flux_kx, label=self.titles[flux])
                ax_log_kx.plot(kx, np.abs(all_flux_kx), ls="--", color=baselogkx.get_color())
                ax_lin_kx.plot(kx, all_flux_kx, label=self.titles[flux])
                baselogky, = ax_log_ky.plot(ky, pos_flux_ky, label=self.titles[flux])
                ax_log_ky.plot(ky, np.abs(all_flux_ky), ls="--", color=baselogky.get_color())
                ax_lin_ky.plot(ky, all_flux_ky, label=self.titles[flux])
                #np.savez("fluxspectra_{}_{}".format(fluxes.species, flux), ky=ky,
                #         fluxspectrum=all_flux_ky)
            ax_log_kx.loglog()
            ax_log_ky.loglog()
            ax_lin_kx.set_xlim(left=0)
            ax_log_kx.set_xlabel(r"$k_x \rho_{ref}$")
            ax_lin_kx.set_xlabel(r"$k_x \rho_{ref}$")
            ax_log_ky.set_xlabel(r"$k_y \rho_{ref}$")
            ax_lin_ky.set_xlabel(r"$k_y \rho_{ref}$")
            for ax in [ax_log_kx, ax_lin_kx, ax_lin_ky, ax_log_ky]:
                # ax.set_ylabel(r"$<|A|^2>$")
                ax.legend()
            ax_log_ky.set_title("{}".format(fluxes.species))
            ax_log_kx.set_title("{}".format(fluxes.species))
            fig_kx.tight_layout()
            fig_ky.tight_layout()
