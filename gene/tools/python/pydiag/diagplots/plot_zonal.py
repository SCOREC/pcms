""" Classes to do plots of the potential and its derivatives """
import numpy as np
import pydiag.utils.averages as av
from pydiag.diagplots.baseplot import plt, Plotting
import pydiag.utils.fourier as fourier
from pydiag.utils import errors as err
from pydiag.utils.geom import Geometry
import pydiag.data.slices
from pydiag.utils.comm import DiagSpace


class ZonalAverage(pydiag.data.slices.MomFieldSlice):
    """ Create a flux-surface averaged ky=0 time series from the field file

    :param common: CommonData object of the run
    :param xrange: Radial range or kx modes to consider
    :param potential: "phi" or "apar"
    """

    def __init__(self, common, rundatafiles, xrange=(None, None), potential="phi", withext=True):
        self.diagspace = DiagSpace(common.spatialgrid, x_fourier=common.x_local, y_fourier=True,
                                   xrange=xrange, yrange=(0, 1), zavg=True)
        super().__init__(common, rundatafiles=rundatafiles, quantity=potential, species=None,
                         diagspace=self.diagspace,
                         modifier_func=lambda dum: dum)  # Use the unmodified data
        self.cm = common
        self.xrange = xrange
        self.phi0_ext = 0  # For the sinusoidal external potential
        self.nk_ext = 1
        self.k_ext = 1
        self.phase_phi_ext = 0
        self.Cxy = 1  # Geometry correction
        # External potential only plausible with electrostatic potential
        self.withext = withext and (potential == "phi")

    def generate_timeseries(self):
        """ Wrapper for fieldlib: Read the ky=0 component of the potential"""
        if not self.cm.y_local:
            raise NotImplementedError('y global is not supported!')
        super().generate_timeseries()
        try:
            phi0_ext = self.cm.pnt.phi0_ext
            nk_ext = self.cm.pnt.kxind_phi_ext
            phase_phi_ext = self.cm.pnt.phase_phi_ext
            k_ext = nk_ext/self.cm.pnt.lx*2*np.pi
        except AttributeError:
            phi0_ext = 0
            nk_ext = 1
            phase_phi_ext = 0
            k_ext = 1
        geom = Geometry(self.cm)
        try:
            self.Cxy = geom.Cxy[self.diagspace.xslice]
        except TypeError:
            self.Cxy = geom.Cxy
        except AttributeError:
            self.Cxy = 1
        if self.withext:
            if self.diagspace.x_fourier:
                xr_min, xr_max = self.xrange[0], self.xrange[1]
                if not self.xrange[0]:
                    xr_min = 0
                if not self.xrange[1]:
                    xr_max = np.inf  # Just a convenient large number
                if xr_min <= nk_ext <= xr_max:
                    nk_ext_subindex = nk_ext - xr_min
                    self.dataarray[:, -nk_ext_subindex] += -0.5j*phi0_ext*np.exp(1j*phase_phi_ext)
                    self.dataarray[:, +nk_ext_subindex] += 0.5j*phi0_ext*np.exp(-1j*phase_phi_ext)
            else:
                xgrid = self.cm.spatialgrid.x[self.diagspace.xslice]
                self.dataarray += phi0_ext*np.sin(k_ext*xgrid + phase_phi_ext)[:, np.newaxis]
        self.dataarray = np.squeeze(self.dataarray)

    def derivative(self, order=1):
        base = self.dataarray[:, self.diagspace.xslice]
        if self.diagspace.x_fourier:
            deriv = (1j*self.cm.spatialgrid.kx_grid_slice(self.diagspace.xslice)) ** order*base
        else:
            xs = self.cm.spatialgrid.x_a[self.diagspace.xslice]
            deriv = np.gradient(base, xs/self.cm.pnt.rhostar, axis=1)
            for o in range(order - 1):
                deriv = np.gradient(deriv, xs/self.cm.pnt.rhostar, axis=1)
        return deriv


class PlotZonalxt(Plotting):
    """ Class to plot the fields from GENE and their derivatives
    So far: electrostatic potential, xt contour

    Only considers the z-averaged ky=0 mode, so far!
    :param commonlist: CommonData objects of the runs to plot
    :param fieldavlist: ZonalAverage objects of the runs to plot
    """

    def __init__(self, commonlist, fieldavlist):
        super().__init__()
        self.fieldavlist = fieldavlist
        self.cm = commonlist
        self.xs = 0

    def plot_phi_xt(self, withext=True, poutput=False):
        """ Plot a color map of phi in x-t space

        :param withext: Include the sinusoidal external potential
        :param poutput: Write plots into pdf files as well
        """
        for cm, fieldav in zip(self.cm, self.fieldavlist):
            fieldav.generate_timeseries()
            if not np.any(fieldav.dataarray):
                print("No field data present. Aborting plot")
                continue
            # Consider x-local and x-global here
            if cm.x_local:
                phi_ky0_kx = fieldav.dataarray
                wexb_t_kx = fieldav.derivative(order=2)/fieldav.Cxy
                self.xs = cm.spatialgrid.x
                phi_ky0 = fourier.kx_to_x(phi_ky0_kx, cm.pnt.nx0, axis=-1)
                wexb_t = fourier.kx_to_x(wexb_t_kx, cm.pnt.nx0, axis=-1)
            else:
                phi_ky0 = np.real_if_close(fieldav.dataarray)[:, fieldav.diagspace.xslice]
                self.xs = cm.spatialgrid.x_a[fieldav.diagspace.xslice]
                wexb_t = -np.real_if_close(-fieldav.derivative(order=2)/fieldav.Cxy)

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            absval = np.max(np.abs(phi_ky0))
            ax1.set_xlim((np.min(self.xs), np.max(self.xs)))
            ax1.set_ylim((np.min(fieldav.timearray), np.max(fieldav.timearray)))
            cm1 = ax1.pcolormesh(self.xs, fieldav.timearray, phi_ky0, shading="gouraud",
                                 cmap=self.cmap_bidirect, vmax=absval, vmin=-absval,
                                 rasterized=True)
            if cm.x_local:
                ax1.set_xlabel(r"$x/\rho_s$")
            else:
                ax1.set_xlabel(r"$x/a$")
            ax1.set_ylabel(r"$t/(a/c_s)$")
            fig1.colorbar(cm1)
            ax1.set_title(r"$\langle\phi\rangle$", y=1.03)
            fig1.tight_layout()

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            absval = np.max(np.abs(wexb_t))
            ax2.set_xlim([np.min(self.xs), np.max(self.xs)])
            ax2.set_ylim([np.min(fieldav.timearray), np.max(fieldav.timearray)])
            cm2 = ax2.pcolormesh(self.xs, fieldav.timearray, wexb_t, shading="gouraud",
                                 cmap=self.cmap_bidirect, vmax=absval, vmin=-absval,
                                 rasterized=True)
            ax2.set_xlabel(r"$x/\rho_s$")
            ax2.set_ylabel(r"$t/(a/c_s)$")
            fig2.colorbar(cm2)
            ax2.set_title(r"$\langle\omega_{E\times B}\rangle$", y=1.03)
            fig2.tight_layout()
            if poutput == 1:
                fig1.savefig('phi_xt{}.pdf'.format(fieldav.cm.fileextension))
                fig2.savefig('wexb_xt{}.pdf'.format(fieldav.cm.fileextension))


class PlotZonalAverage(Plotting):
    """ Diagnostic to investigate zonal flows (the ky=0 mode of the potential)

    :param commonlist: CommonData objects of the runs to plot
    :param fieldavlist: ZonalAverage objects of the runs to plot
    """

    def __init__(self, commonlist, fieldavlist):
        super().__init__()
        self.fieldavlist = fieldavlist
        self.xs = 0
        self.kx = 0
        self.cm = commonlist

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
        for cm, fieldav in zip(self.cm, self.fieldavlist):
            fieldav.generate_timeseries()
            if not np.any(fieldav.dataarray):
                print("No field data present. Aborting calculation and plot")
                continue
            # Consider x-local and x-global here
            if cm.x_local:
                self.kx = cm.spatialgrid.kx_grid_slice(fieldav.diagspace.xslice)
                self.rms_omega_local(fieldav, axlist)
            else:
                self.xs = cm.spatialgrid.x[fieldav.diagspace.xslice]
                self.xs = self.xs*cm.pnt.rhostar + cm.pnt.x0
                self.rms_omega_global(fieldav, axlist)
        for fig in figlist:
            fig.tight_layout()
        if poutput:
            figlist[0].savefig('phi_zonal{}.pdf'.format(self.fieldavlist[-1].cm.fileextension))
            figlist[1].savefig('erad_zonal{}.pdf'.format(self.fieldavlist[-1].cm.fileextension))
            figlist[2].savefig('wexb_zonal{}.pdf'.format(self.fieldavlist[-1].cm.fileextension))

    def rms_omega_local(self, fieldav, axlist):
        """ Local version for the ExB shear calculation

        Generates time averaged kx spectra for <phi>, E_r and omega_ExB
        and prints the RMS for omega (order: first x RMS then time average)
        The plots are logarithmic and ignore the kx=0 mode
        (would be NC, but we cannot calculate phi for it anyway),
        as well as the negative kx - which are Hermitian to kx>0 for ky=0 anyway

        :param fieldav: The ZonalAverage object to be plotted
        :param axlist: The list of axes objects to plot into (len(axlist) >= 3)
        """
        kx = self.kx
        phi_ky0 = fieldav.dataarray
        erad_t = -fieldav.derivative(order=1)
        omega_t = -fieldav.derivative(order=2)/fieldav.Cxy
        times = fieldav.timearray
        rms_omega_t = np.sqrt(np.sum(np.square(np.abs(omega_t)), axis=1))
        err_rms_omega = err.windowerr(rms_omega_t, times)[0]
        av_rms_omega = av.mytrapz(rms_omega_t, times)
        print("Rms ExB shear (time average): ", av_rms_omega, " +- ", err_rms_omega)
        phi_tav = av.mytrapz(np.fft.ifftshift(np.abs(phi_ky0)), times)
        erad_tav = av.mytrapz(np.fft.ifftshift(np.abs(erad_t)), times)
        omega_tav = av.mytrapz(np.fft.ifftshift(np.abs(omega_t)), times)

        ax1 = axlist[0]
        ax1.set_xscale('log')
        ax1.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax1.set_yscale('log')
        ax1.set_xlabel(r'$k_x \rho_s$')
        ax1.set_ylabel(r"$|\langle \phi \rangle|$")
        ax1.plot(np.fft.ifftshift(kx), phi_tav)
        ax1.set_ylim(ymin=self._minlogaxis(phi_tav), ymax=self._maxlogaxis(phi_tav))

        ax2 = axlist[1]
        ax2.set_xscale('log')
        ax2.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$k_x \rho_s$')
        ax2.set_ylabel(r"$|\langle \mathbf{E_x} \rangle|$")
        ax2.plot(np.fft.ifftshift(kx), erad_tav)
        ax2.set_ylim(ymin=self._minlogaxis(erad_tav), ymax=self._maxlogaxis(erad_tav))

        ax3 = axlist[2]
        ax3.set_xscale('log')
        ax3.set_xlim(xmin=self._minlogaxis(np.abs(kx)), xmax=self._maxlogaxis(np.abs(kx)))
        ax3.set_yscale('log')
        ax3.set_xlabel(r'$k_x \rho_s$')
        ax3.set_ylabel(r"$|\langle \omega_{\mathbf{E}\times\mathbf{B}} \rangle|$")
        ax3.plot(np.fft.ifftshift(kx), omega_tav)
        ax3.set_ylim(ymin=self._minlogaxis(omega_tav), ymax=self._maxlogaxis(omega_tav))

    def rms_omega_global(self, fieldav, axlist):
        """ Global version for the ExB shear calculation

        Generates time averaged x profiles for <phi>, E_r and omega_ExB
        and prints the RMS and median of absolutes
        for omega (order: first x RMS then time average)

        :param fieldav: The ZonalAverage object to be plotted
        :param axlist: The list of axes objects to plot into (len(axlist) >= 3)
        """
        xs = self.xs
        phi_ky0_x = np.real_if_close(fieldav.dataarray)
        erad_t = -np.real_if_close(fieldav.derivative(order=1))
        wexb_t = -np.real_if_close(fieldav.derivative(order=2)/fieldav.Cxy)
        fullphi_tav = av.mytrapz(phi_ky0_x, fieldav.timearray)
        erad_tav = av.mytrapz(erad_t, fieldav.timearray)
        wexb_tav = av.mytrapz(wexb_t, fieldav.timearray)

        rms_wexb_t = np.sqrt(np.mean(wexb_t ** 2, axis=1))
        err_rms_wexb = err.windowerr(rms_wexb_t, np.array(fieldav.timearray))[0]
        av_rms_wexb = av.mytrapz(rms_wexb_t, fieldav.timearray)
        print("Rms ExB shear (time average): ", av_rms_wexb, " +- ", err_rms_wexb)
        medabs_wexb_t = np.median(np.abs(wexb_t), axis=1)
        av_medabs_wexb = av.mytrapz(medabs_wexb_t, fieldav.timearray)
        print("Median absolute ExB shear (time average): ", av_medabs_wexb)

        ax1 = axlist[0]
        ax1.set_xlabel(r'$x/a$')
        ax1.set_ylabel(r"$\langle \phi \rangle$")
        ax1.plot(xs, fullphi_tav)

        ax2 = axlist[1]
        ax2.set_xlabel(r'$x/a$')
        ax2.set_ylabel(r"$\langle \mathbf{E_x} \rangle$")
        ax2.plot(xs, erad_tav)

        ax3 = axlist[2]
        ax3.set_xlabel(r'$x/a$')
        ax3.set_ylabel(r"$\langle \omega_{\mathbf{E}\times\mathbf{B}} \rangle$")
        ax3.plot(xs, wexb_tav)
        ax3.axhline(color="Black", ls="--", label='_nolegend_')
