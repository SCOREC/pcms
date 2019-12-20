""" Classes to do plots of the potential and its derivatives """
import numpy as np

import pydiag.utils.averages as av
from pydiag.diagplots.baseplot import plt, Plotting
import pydiag.utils.fourier as fourier
from pydiag.utils import errors as err


class PlotField(Plotting):
    """ Class to plot the fields from GENE and their derivatives
    So far: electrostatic potential, xt contour

    Only considers the z-averaged ky=0 mode, so far!
    :param commonlist: CommonData objects of the runs to plot
    :param phiflist: FieldFile objects of the runs to plot
    """

    def __init__(self, commonlist, phiflist):
        super().__init__()
        self.phiflist = phiflist
        self.cm = commonlist
        self.xs = 0

    def plot_phi_xt(self, withext=True, poutput=False):
        """ Plot a color map of phi in x-t space

        :param withext: Include the sinusoidal external potential
        :param poutput: Write plots into pdf files as well
        """
        for cm, phif in zip(self.cm, self.phiflist):
            if not np.any(phif.dataarray):
                print("No field data present. Aborting plot")
                continue
            # Consider x-local and x-global here
            if cm.x_local:
                phi_ky0_kx = np.copy(phif.dataarray)
                if withext:
                    phi_ky0_kx[:, -phif.nk_ext] += -0.5j*phif.phi0_ext*np.exp(1j*phif.phase_phi_ext)
                    phi_ky0_kx[:, +phif.nk_ext] += 0.5j*phif.phi0_ext*np.exp(-1j*phif.phase_phi_ext)
                wexb_t_kx = - cm.spatialgrid.kx_fftorder ** 2/phif.Cxy*phi_ky0_kx
                self.xs = cm.spatialgrid.x
                phi_ky0 = fourier.kx_to_x(phi_ky0_kx, cm.pnt.nx0, axis=-1)
                wexb_t = fourier.kx_to_x(wexb_t_kx, cm.pnt.nx0, axis=-1)
            else:
                phi_ky0 = np.real_if_close(phif.dataarray)
                self.xs = cm.spatialgrid.x_a
                erad_t = np.empty(phi_ky0.shape)
                wexb_t = np.empty(phi_ky0.shape)
                for tind in range(len(phif.timearray)):
                    erad_t[tind, :] = -np.gradient(phi_ky0[tind, :],
                                                   (self.xs[1]-self.xs[0])/cm.pnt.rhostar)
                    wexb_t[tind, :] = np.gradient(erad_t[tind, :]/phif.Cxy,
                                                  (self.xs[1]-self.xs[0])/cm.pnt.rhostar)

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            absval = np.max(np.abs(phi_ky0))
            ax1.set_xlim([np.min(self.xs), np.max(self.xs)])
            ax1.set_ylim([np.min(phif.timearray), np.max(phif.timearray)])
            cm1 = ax1.pcolormesh(self.xs, phif.timearray, phi_ky0, shading="gouraud",
                                 cmap=self.cmap_bidirect, vmax=absval, vmin=-absval, rasterized=True)
            ax1.set_xlabel(r"$x/\rho_s$", fontsize=self.xyfs)
            ax1.set_ylabel(r"$t/(a/c_s)$", fontsize=self.xyfs)
            fig1.colorbar(cm1)
            ax1.set_title(r"$\langle\phi\rangle$", fontsize=self.titlefs, y=1.03)
            fig1.tight_layout()

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            absval = np.max(np.abs(wexb_t))
            ax2.set_xlim([np.min(self.xs), np.max(self.xs)])
            ax2.set_ylim([np.min(phif.timearray), np.max(phif.timearray)])
            cm2 = ax2.pcolormesh(self.xs, phif.timearray, wexb_t, shading="gouraud",
                                 cmap=self.cmap_bidirect, vmax=absval, vmin=-absval, rasterized=True)
            ax2.set_xlabel(r"$x/\rho_s$", fontsize=self.xyfs)
            ax2.set_ylabel(r"$t/(a/c_s)$", fontsize=self.xyfs)
            fig2.colorbar(cm2)
            ax2.set_title(r"$\langle\omega_{E\times B}\rangle$", fontsize=self.titlefs, y=1.03)
            fig2.tight_layout()
            if poutput == 1:
                fig1.frameon = False
                fig2.frameon = False
                fig1.savefig('phi_xt{}.pdf'.format(phif.cm.fext))
                fig2.savefig('wexb_xt{}.pdf'.format(phif.cm.fext))

    @classmethod
    def show(cls):
        plt.show()


class Zonal(Plotting):
    """ Diagnostic to investigate zonal flows (the ky=0 mode of the potential)

    :param commonlist: CommonData objects of the runs to plot
    :param phiflist: FieldFile objects of the runs to plot
    """

    def __init__(self, commonlist, phiflist):
        super().__init__()
        self.phiflist = phiflist
        self.xs = 0
        self.kx = 0
        self.cm = commonlist

    @staticmethod
    def _minlogaxis(arr):
        """ Calculate lower limit for plotting: lower Integer * 10^minimumexponent

        :param arr: Array to base the limit on
        """
        minim = arr[arr != 0].min()  # ignore 0
        min_exp = np.floor(np.log10(minim))
        min_mant = np.floor(minim*10**(-min_exp))
        floormin = min_mant*10**(min_exp)
        return floormin

    @staticmethod
    def _maxlogaxis(arr):
        """ Calculate upper limit for plotting: higher Integer * 10^maximumexponent

        :param arr: Array to base the limit on
        """
        maxim = arr[arr != 0].max()
        max_exp = np.floor(np.log10(maxim))
        max_mant = np.ceil(maxim*10**(-max_exp))
        ceilmax = max_mant*10**(max_exp)
        return ceilmax

    def createfigs(self, withext=False, poutput=False):
        """ Plot the time averaged root mean square of the ExB shear

        This routine separates between the local and global shear calculation and
        handles the general frameworks for the plots
        (all lines in one plot or multiple windows)
        :param withext: Include the external sinusoidal potential
        :param poutput: Generate pdf files of the plots
        """
        figlist = [plt.figure() for i in range(3)]
        axlist = [fig.add_subplot(111) for fig in figlist]
        for cm, phif in zip(self.cm, self.phiflist):
            if not np.any(phif.dataarray):
                print("No field data present. Aborting calculation and plot")
                continue
            # Consider x-local and x-global here
            if cm.x_local:
                self.kx = 2*np.pi * np.fft.fftfreq(cm.pnt.nx0, d=cm.pnt.lx/cm.pnt.nx0)
                self.rms_omega_local(phif, axlist, withext)
            else:
                self.xs = np.array([-cm.pnt.lx/2.+cm.pnt.lx/(cm.pnt.nx0-1)*i
                                    for i in range(0, cm.pnt.nx0)])
                self.xs = self.xs * cm.pnt.rhostar + cm.pnt.x0
                self.rms_omega_global(phif, axlist, withext)
        for fig in figlist:
            fig.tight_layout()
        if poutput:
            figlist[0].frameon = False
            figlist[0].savefig('phi_zonal{}.pdf'.format(self.phiflist[-1].cm.fext))
            figlist[1].frameon = False
            figlist[1].savefig('erad_zonal{}.pdf'.format(self.phiflist[-1].cm.fext))
            figlist[2].frameon = False
            figlist[2].savefig('wexb_zonal{}.pdf'.format(self.phiflist[-1].cm.fext))

    def rms_omega_local(self, phif, axlist, withext=False):
        """ Local version for the ExB shear calculation

        Generates time averaged kx spectra for <phi>, E_r and omega_ExB
        and prints the RMS for omega (order: first x RMS then time average)
        The plots are logarithmic and ignore the kx=0 mode
        (would be NC, but we cannot calculate phi for it anyway),
        as well as the negative kx - which are Hermitian to kx>0 for ky=0 anyway

        :param phif: The FieldData object to be plotted
        :param axlist: The list of axes objects to plot into (len(axlist) >= 3)
        :param withext: Include the external sinusoidal potential
        """
        phi_ky0 = np.copy(phif.dataarray)
        kx = self.kx
        if withext:
            phi_ky0[:, -phif.nk_ext] += -0.5j*phif.phi0_ext*np.exp(1j*phif.phase_phi_ext)
            phi_ky0[:, +phif.nk_ext] += 0.5j*phif.phi0_ext*np.exp(-1j*phif.phase_phi_ext)
        erad_t = - kx * phi_ky0
        omega_t = - kx**2/phif.Cxy * phi_ky0
        times = phif.timearray
        rms_omega_t = np.sqrt(np.sum(np.square(np.abs(omega_t)), axis=1))
        err_rms_omega = err.windowerr(rms_omega_t, times)[0]
        av_rms_omega = av.mytrapz(rms_omega_t, times)
        print("Rms ExB shear (time average): ", av_rms_omega,
              " +- ", err_rms_omega)
        phi_tav = av.mytrapz(np.fft.ifftshift(np.abs(phi_ky0)), times)
        erad_tav = av.mytrapz(np.fft.ifftshift(np.abs(erad_t)), times)
        omega_tav = av.mytrapz(np.fft.ifftshift(np.abs(omega_t)), times)

        ax1 = axlist[0]
        ax1.set_xscale('log')
        ax1.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax1.set_yscale('log')
        ax1.set_xlabel(r'$k_x \rho_s$', fontsize=self.xyfs)
        ax1.set_ylabel(r"$|\langle \phi \rangle|$", fontsize=self.xyfs)
        ax1.plot(np.fft.ifftshift(kx), phi_tav)
        ax1.set_ylim(ymin=self._minlogaxis(phi_tav), ymax=self._maxlogaxis(phi_tav))

        ax2 = axlist[1]
        ax2.set_xscale('log')
        ax2.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$k_x \rho_s$', fontsize=self.xyfs)
        ax2.set_ylabel(r"$|\langle \mathbf{E_x} \rangle|$", fontsize=self.xyfs)
        ax2.plot(np.fft.ifftshift(kx), erad_tav)
        ax2.set_ylim(ymin=self._minlogaxis(erad_tav), ymax=self._maxlogaxis(erad_tav))

        ax3 = axlist[2]
        ax3.set_xscale('log')
        ax3.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax3.set_yscale('log')
        ax3.set_xlabel(r'$k_x \rho_s$', fontsize=self.xyfs)
        ax3.set_ylabel(r"$|\langle \omega_{\mathbf{E}\times\mathbf{B}} \rangle|$",
                       fontsize=self.xyfs)
        ax3.plot(np.fft.ifftshift(kx), omega_tav)
        ax3.set_ylim(ymin=self._minlogaxis(omega_tav), ymax=self._maxlogaxis(omega_tav))

    def rms_omega_global(self, phif, axlist, withext=False):
        """ Global version for the ExB shear calculation

        Generates time averaged x profiles for <phi>, E_r and omega_ExB
        and prints the RMS and median of absolutes
        for omega (order: first x RMS then time average)

        :param phif: The fielddata object to be plotted
        :param axlist: The list of axes objects to plot into (len(axlist) >= 3)
        :param withext: Include the external sinusoidal potential
        """
        xs = self.xs
        phi_ky0_x = np.real_if_close(phif.dataarray)
        if withext:  # This creates a full copy of the phi array
            fullphi = phi_ky0_x + phif.phi0_ext*np.sin(phif.k_ext*xs + phif.phase_phi_ext)
        else:
            fullphi = phi_ky0_x
        erad_t = np.empty(fullphi.shape)
        wexb_t = np.empty(fullphi.shape)
        for tind in range(len(phif.timearray)):
            erad_t[tind, :] = -np.gradient(fullphi[tind, :], (xs[1]-xs[0])/phif.cm.pnt.rhostar)
            wexb_t[tind, :] = np.gradient(erad_t[tind, :]/phif.Cxy, (xs[1]-xs[0])/phif.cm.pnt.rhostar)
        fullphi_tav = av.mytrapz(fullphi, phif.timearray)
        erad_tav = av.mytrapz(erad_t, phif.timearray)
        wexb_tav = av.mytrapz(wexb_t, phif.timearray)

        rms_wexb_t = np.sqrt(np.mean(wexb_t**2, axis=1))
        err_rms_wexb = err.windowerr(rms_wexb_t, np.array(phif.timearray))[0]
        av_rms_wexb = av.mytrapz(rms_wexb_t, phif.timearray)
        print("Rms ExB shear (time average): ", av_rms_wexb, " +- ", err_rms_wexb)
        medabs_wexb_t = np.median(np.abs(wexb_t), axis=1)
        av_medabs_wexb = av.mytrapz(medabs_wexb_t, phif.timearray)
        print("Median absolute ExB shear (time average): ", av_medabs_wexb)

        ax1 = axlist[0]
        ax1.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax1.set_ylabel(r"$\langle \phi \rangle$", fontsize=self.xyfs)
        ax1.plot(xs, fullphi_tav)

        ax2 = axlist[1]
        ax2.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax2.set_ylabel(r"$\langle \mathbf{E_x} \rangle$", fontsize=self.xyfs)
        ax2.plot(xs, erad_tav)

        ax3 = axlist[2]
        ax3.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax3.set_ylabel(r"$\langle \omega_{\mathbf{E}\times\mathbf{B}} \rangle$",
                       fontsize=self.xyfs)
        ax3.plot(xs, wexb_tav)
        ax3.axhline(color="Black", ls="--", label='_nolegend_')

    @classmethod
    def show(cls):
        plt.show()
