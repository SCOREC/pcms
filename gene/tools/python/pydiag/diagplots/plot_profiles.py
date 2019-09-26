"""PlotProfiles.py: Module to do the plotting of profile diagnostic outputput """

import numpy as np
import pydiag.utils.averages as av
from pydiag.utils import errors as err, nc_predictions
from pydiag.diagplots.baseplot import plt, Plotting


def printaverageflux(prdlist, xmin, xmax):
    """ Print the (time averaged) mean heat flux for a radial domain

    :param prdlist: List of ProfileData objects to process
    :param xmin: Left boundary in x/L_ref
    :param xmax: Right boundary in x/L_ref
    """
    for prd in prdlist:
        for prs in prd.prspec:
            meanQnc = av.mytrapz(prs.dataarray["Qnc"], prs.timearray)
            errorQnc = err.windowerr(prs.dataarray["Qnc"], prs.timearray)[0]
            errorQt = err.windowerr(prs.dataarray["Qturb"], prs.timearray)[0]
            meanQt = av.mytrapz(prs.dataarray["Qturb"], prs.timearray)
            meanomn = av.mytrapz(prs.dataarray["omns"], prs.timearray)
            meanomt = av.mytrapz(prs.dataarray["omts"], prs.timearray)
            # Calculate the average fluxes at x=xmin-xmax and print them
            postoav = np.where(np.logical_and(prd.xs > xmin, prd.xs < xmax))
            meanQt_center = np.mean(meanQt[postoav])
            errorQt_center = np.mean(errorQt[postoav])
            meanomn_center = np.mean(meanomn[postoav])
            meanomt_center = np.mean(meanomt[postoav])
            print("Mean turbulent flux at {:2}-{:2}: {:5} +- {:5}".format(xmin, xmax, meanQt_center,
                                                                          errorQt_center))
            print("Mean density gradient at {:2}-{:2}: {:5}".format(xmin, xmax, meanomn_center))
            print("Mean temp gradient at {:2}-{:2}: {:5}".format(xmin, xmax, meanomt_center))
            meanQnc_center = np.mean(meanQnc[postoav])
            errorQnc_center = np.mean(errorQnc[postoav])
            print("Mean NC flux at {:2}-{:2}: {:5} +- {:5}".format(xmin, xmax, meanQnc_center,
                                                                   errorQnc_center))


def calc_fb_param(prd, fieldd):
    """Calculate the k in the radial force balance equation

    Uses time averaged quantities.
    :param prd: ProfileData object for n, T and the parameters file
    :param fieldd: FieldData object for E_r
    :returns: numpy array of k for each radial position
    """
    phi_ky0 = np.atleast_3d(np.real_if_close(fieldd.dataarray))
    erad_t = np.empty(phi_ky0.shape)
    for tind in range(len(fieldd.timearray)):
        erad_t[tind, :, 0] = -np.gradient(phi_ky0[tind, :, 0],
                                       (prd.xs[1] - prd.xs[0])/fieldd.cm.pnt.rhostar)
    erad = av.mytrapz(erad_t[..., 0], fieldd.timearray)
    fb_param = []
    epsilon = np.array(prd.xs*prd.cm.pnt.minor_r/prd.cm.pnt.major_R)
    btoverbp = prd.q/epsilon*np.sqrt(1 - epsilon ** 2)
    for prs in prd.prspec:
        times = prs.timearray
        jb_av = av.mytrapz(prs.dataarray["jbs"], times)
        ns_av = av.mytrapz(prs.dataarray["ns"], times)
        Ts_av = av.mytrapz(prs.dataarray["Ts"], times)
        omns_av = av.mytrapz(prs.dataarray["omns"], times)
        omts_av = av.mytrapz(prs.dataarray["omts"], times)
        fb_param.append(-(jb_av/(btoverbp*ns_av*Ts_av) - (omns_av + omts_av + erad/Ts_av))/omts_av)
    return fb_param


class PlotProfdata(Plotting):
    """ Plot time averaged profiles

    Generates plots of temperature and density, particle, heat and momentum flux,
    bootstrap current, radial force balance parameter, collisionality
    :param profdata: ProfileData object list to plot
    :param fielddata: FieldData object list to use in force balance plots
    :param plotset: Set which controls the plots
    """

    def __init__(self, profdata, fielddata, plotset):
        super().__init__()
        # Data structures used in the different plots
        self.prdlist = profdata
        self.fieldlist = fielddata
        self.plotset = plotset
        for prd in self.prdlist:
            if prd.isdatapresent is False:
                print("No profile data present in object ", prd.fext)
                self.plotset.clear()
        # Standard limits for the x axis are determined over all input runs
        self.xmin = min([prd.xs[1] for prd in self.prdlist])
        self.xmax = max([prd.xs[-1] for prd in self.prdlist])

    def output_avprof(self, avproffile_prefix):
        """ Generate time averaged density and temperature profiles
        which can be used as input for GENE

        :param avproffile_prefix: Output filename, species name is added to it
        """
        print("Writing averaged profile to file")
        for prd in self.prdlist:
            for n, prs in enumerate(prd.prspec):
                av_T = av.mytrapz(prs.Ts, prs.timearray)
                av_n = av.mytrapz(prs.ns, prs.timearray)
                avprofilename = (avproffile_prefix + "_" + prd.specname[n])
                with open(avprofilename, 'w') as avproffile:
                    avproffile.write(
                            "#   x/a             x/a(dummy)       T/Tref          n/nref\n")
                    avproffile.write("# Averaged profile generated by profile_data.py:"
                                     "t = {:4.1f} - {:4.1f}\n".format(prd.cm.pnt.starttime,
                                                                      prd.cm.pnt.endtime))
                    for i in range(0, prd.pnt.nx0):
                        avproffile.write(
                                '{:E}'.format(prd.xs[i]) + "\t" + '{:E}'.format(prd.xs[i]) + "\t")
                        avproffile.write(
                                '{:E}'.format(av_T[i]) + "\t" + '{:E}'.format(av_n[i]) + "\n")

    def createfigs(self, poutput):
        """ Generate the figures and axes objects to plot into for averaged profiles

        Calls the individual plotting routines
        Plots individual figures for each quantity and species, but puts their data from different
        runs into the same figures
        :param poutput: Switch to generate pdf file output
        """
        # Number of plotted variables
        numvars = 5
        # Number of plot windows
        nspecs = self.prdlist[0].cm.pnt.n_spec
        numfigs = (numvars - 1)*nspecs + 2
        figlist = [plt.figure() for _ in range(numfigs)]
        # One plot window for each species and quantity, one axis object in each
        axlist = [fig.add_subplot(111) for fig in figlist]
        # One legend per subplot
        leglist = [[] for _ in axlist]
        for iprd, prd in enumerate(self.prdlist):
            # Use the list of axes so that the different ProfileData objects end up in
            # the same figure and axis object (one per species)
            if 'plotT' in self.plotset:
                self.plot_tempdens(axlist[0:nspecs], prd)
            if 'plotQ' in self.plotset or 'plotG' in self.plotset or 'plotP' in self.plotset:
                self.plot_radflux(axlist[nspecs:2*nspecs], prd, "Qturb")
                self.plot_radflux(axlist[nspecs:2*nspecs], prd, "Qnc")
                # NC predictions only for first species (ions)
                self.plot_Qncpred(axlist[nspecs], prd)
                self.plot_radflux_SI(axlist[2*nspecs:3*nspecs], prd, "Qturb")
                self.plot_radflux_SI(axlist[2*nspecs:3*nspecs], prd, "Qnc")
            if 'plotj' in self.plotset:
                # Exception: Put all species in one axis
                self.plot_jbs([axlist[3*nspecs] for _ in range(nspecs)], prd, self.fieldlist[iprd])
                self.plot_fb_param([axlist[3*nspecs + 1] for _ in range(nspecs)], prd,
                                   self.fieldlist[iprd])
            axcount = 3*nspecs + 2
            if 'plotnu' in self.plotset:
                self.plot_nustar(axlist[axcount:axcount + nspecs], prd)
            self.shade_buffer_region(axlist, prd)
        for ax in axlist:
            ax.legend()
        if poutput:
            figlist[0].savefig("profiles{}.pdf".format(self.prdlist[0].fext))
            figlist[1].savefig("radfluxes{}.pdf".format(self.prdlist[0].fext))
            figlist[2].savefig("radareaflux{}.pdf".format(self.prdlist[0].fext))
            figlist[3].savefig("j_botstrap{}.pdf".format(self.prdlist[0].fext))
            figlist[4].savefig("nustar{}.pdf".format(self.prdlist[0].fext))

    @staticmethod
    def shade_buffer_region(axlist, prd):
        ubuff = 0
        lbuff = 0
        try:
            if prd.pnt.lcoef_krook:
                lbuff = prd.pnt.l_buffer_size
        except AttributeError:  # No buffer zone set
            pass
        try:
            if prd.pnt.ucoef_krook:
                ubuff = prd.pnt.u_buffer_size
        except AttributeError:
            pass
        for ax in axlist:
            if lbuff > 0:
                ax.axvspan(prd.xs[0], prd.xs[0] + lbuff, alpha=0.4, color="Black")
            if ubuff > 0:
                ax.axvspan(prd.xs[-1] - ubuff, prd.xs[-1], alpha=0.4, color="Black")

    def plot_tempdens(self, axlist, prd):
        """Plot the temperature and density profiles

        Profiles are temporally averaged by the trapezoid rule
        if the lists contain more than one time block.
        The initial profiles are given as reference.

        :param axlist: List of axes objects to plot into (one per species)
        :param prd: ProfileData object containing the data
        """
        for iprs, prs in enumerate(prd.prspec):
            baseomn, = axlist[iprs].plot(prd.xs, av.mytrapz(prs.dataarray["omns"], prs.timearray),
                                         label=r'$a/L_n$')
            baseomt, = axlist[iprs].plot(prd.xs, av.mytrapz(prs.dataarray["omts"], prs.timearray),
                                         label=r'$a/L_T$')
            basen, = axlist[iprs].plot(prd.xs, av.mytrapz(prs.dataarray["ns"], prs.timearray),
                                       label=r'$n$')
            baset, = axlist[iprs].plot(prd.xs, av.mytrapz(prs.dataarray["Ts"], prs.timearray),
                                       label=r'$T$')
            axlist[iprs].plot(prd.xs, prd.omn0s[:, 0], '--', color=baseomn.get_color())
            axlist[iprs].plot(prd.xs, prd.omt0s[:, 0], '--', color=baseomt.get_color())
            axlist[iprs].plot(prd.xs, prd.T0s[:, 0]*prd.cm.pars["temp" + str(iprs + 1)], '--',
                              color=baset.get_color())
            axlist[iprs].plot(prd.xs, prd.n0s[:, 0]*prd.cm.pars["dens" + str(iprs + 1)], '--',
                              color=basen.get_color())
            axlist[iprs].set_xlabel(r'$x/a$')
            axlist[iprs].set_xlim(self.xmin, self.xmax)
            axlist[iprs].set_title(prd.cm.specnames[iprs])

    def plot_radflux(self, axlist, prd, quantity, linestyle="-"):
        """Plot the averaged fluxes

         Plots the time averaged particle, heat and momentum fluxes.

         :param axlist: List of axes objects (one per each species)
         :param prd: ProfileData object containing the data
         :param quantity: Name of the flux to plot
         :param linestyle: See matplotlib documentation for possible forms
         """
        linelabel = {"Gammaturb": r"$\Gamma_{{\mathrm{{turb}}}}$",
                     "Gammanc": r"$\Gamma_{{\mathrm{{NC}}}}$", "Qturb": r"$Q_{{\mathrm{{turb}}}}$",
                     "Qnc": r"$Q_{{\mathrm{{NC}}}}$", "Piturb": r"$\Pi_{{\mathrm{{turb}}}}$",
                     "Pinc": r"$\Pi_{{\mathrm{{NC}}}}$"}
        for n, prs in enumerate(prd.prspec):
            mean_quant = av.mytrapz(prs.dataarray[quantity], prs.timearray)
            error_quant = err.windowerr(prs.dataarray[quantity], prs.timearray)[0]
            base_line, = axlist[n].plot(prd.xs, mean_quant, linestyle=linestyle,
                                        label=linelabel[quantity] + "{}".format(
                                                prd.cm.specnames[n]))
            axlist[n].fill_between(prd.xs, mean_quant - error_quant, mean_quant + error_quant,
                                   alpha=0.2, facecolor=base_line.get_color())
            axlist[n].set_xlim(self.xmin, self.xmax)
            axlist[n].set_xlabel(r'$x/a$')
            self._set_flux_y_label(axlist[n], quantity)

    def _set_flux_y_label(self, axesobj, quantity, SI=False):
        """ Set the y axis label to the right type of gyro-Bohm or SI unit"""
        if SI:
            gbaxlabel = {"Gammaturb": r"$10^19/\mathrm{m^2s}$", "Gammanc": r"$10^19/\mathrm{m^2s}$",
                         "Qturb": r"$\mathrm{kW/m^2}$", "Qnc": r"$\mathrm{kW/m^2}$",
                         "Piturb": r"$\mathrm{kg/m s^2}$", "Pinc": r"$\mathrm{kg/m s^2}$"}
        else:
            gbaxlabel = {"Gammaturb": r"$\Gamma/\Gamma_{gB}$", "Gammanc": r"$\Gamma/\Gamma_{gB}$",
                         "Qturb": r"$Q/Q_{gB}$", "Qnc": r"$Q/Q_{gB}$", "Piturb": r"$\Pi/\Pi_{gB}$",
                         "Pinc": r"$\Pi/\Pi_{gB}$"}
        currlabel = axesobj.get_ylabel()
        if gbaxlabel[quantity] not in currlabel:
            currlabel += gbaxlabel[quantity]
        axesobj.set_ylabel(currlabel)

    @staticmethod
    def plot_Qncpred(axis, prd, potato=False):
        """ Plot neoclassical predictions based on several (semi-)analytic formulas

        Currently provides Chang-Hinton and potato region models from Lin and Bergmann.
        :param axis: A single axes object (should be the main ions, Chang-Hinton might not apply
        for anything else)
        :param prd:  ProfileData object containing the data
        :param potato: Include potato regime predictions
        """
        meanch = av.mytrapz(prd.Qcharr, prd.prspec[0].timearray)
        axis.plot(prd.xs, meanch, '--', label=r"Chang-Hinton")
        if potato:
            potwidth = nc_predictions.calc_potwidth(prd)
            x0 = np.array(0.25*(prd.xs/potwidth*(1 + np.sqrt(5))) ** 3)
            axis.plot(prd.xs, meanch*(1 - np.exp(-x0)*(x0 + 1)), '-.', label=r'Lin')
            rpind3 = np.where(prd.xs <= 3*potwidth)
            axis.plot(prd.xs[rpind3], meanch[rpind3]*(1 - (1 - x0[rpind3] ** (1.0/3)/3) ** 2),
                      label=r'Bergmann')

    def plot_radflux_SI(self, axlist, prd, quantity, linestyle="-"):
        """Plot the averaged fluxes in SI units

         Plots the time averaged particle, heat and momentum fluxes.

         :param axlist: List of axes objects (one per each species)
         :param prd: ProfileData object containing the data
         :param quantity: Name of the flux to plot
         :param linestyle: See matplotlib documentation for possible forms
         """
        linelabel = {"Gammaturb": r"$\Gamma_{{\mathrm{{turb}}}}$",
                     "Gammanc": r"$\Gamma_{{\mathrm{{NC}}}}$", "Qturb": r"$Q_{{\mathrm{{turb}}}}$",
                     "Qnc": r"$Q_{{\mathrm{{NC}}}}$", "Piturb": r"$\Pi_{{\mathrm{{turb}}}}$",
                     "Pinc": r"$\Pi_{{\mathrm{{NC}}}}$"}
        for n, prs in enumerate(prd.prspec):
            mean_quant = av.mytrapz(prs.dataarray[quantity], prs.timearray)
            mean_quant_SI = mean_quant*self.gyrobohm_SI(prd.cm, quantity)
            error_quant = err.windowerr(prs.dataarray[quantity], prs.timearray)[0]
            error_quant_SI = error_quant*self.gyrobohm_SI(prd.cm, quantity)
            base_line, = axlist[n].plot(prd.xs, mean_quant_SI, linestyle=linestyle,
                                        label=linelabel[quantity] + "{}".format(
                                                prd.cm.specnames[n]))
            axlist[n].fill_between(prd.xs, mean_quant_SI - error_quant_SI,
                                   mean_quant_SI + error_quant_SI, alpha=0.2,
                                   facecolor=base_line.get_color())
            axlist[n].set_xlim(self.xmin, self.xmax)
            axlist[n].set_xlabel(r'$x/a$')
            self._set_flux_y_label(axlist[n], quantity, SI=True)

    @staticmethod
    def renorm_profile(prd, ispec, quantity):
        """ Renormalization factor for a radial profile to local gyro-Bohm units

        i.e. gB units of each radial position instead of the gB at x0 position
        :param prd: ProfileData object to get temperature and density from
        :param ispec: Species index for prof
        :param quantity: Type of prof (e.g. Gammaturb, Qnc)
        """
        renorm = {"Ts": prd.T0s[:, ispec], "ns": prd.n0s[:, ispec], "omts": 1,
                  # The T,n  dependency in omt, omn compensates itself
                  "omns": 1, "Gammanc": prd.T0s[:, ispec] ** (3/2)*prd.n0s[:, ispec],
                  "Gammaturb": prd.T0s[:, ispec] ** (3/2)*prd.n0s[:, ispec],
                  "Qturb": prd.T0s[:, ispec] ** (5/2)*prd.n0s[:, ispec],
                  "Qnc": prd.T0s[:, ispec] ** (5/2)*prd.n0s[:, ispec],
                  "Pinc": prd.T0s[:, ispec] ** 2*prd.n0s[:, ispec],
                  "Piturb": prd.T0s[:, ispec] ** 2*prd.n0s[:, ispec],
                  "jbs": prd.T0s[:, ispec]*prd.n0s[:, ispec]}
        return renorm[quantity]

    def plot_jbs(self, axlist, prd, fieldd):
        """Plot the species-wise bootstrap current

        Plots time averaged < u_par B >_FS for every species.
        The bootstrap current would be their sum
        :param axlist: List of axes objects (one per each species)
        :param prd: ProfileData object containing the data
        :param fieldd: FieldData object for radial force balance
        """
        for n, prs in enumerate(prd.prspec):
            axlist[n].plot(prd.xs, av.mytrapz(prs.dataarray["jbs"], prs.timearray),
                           label=r"$j_{BS}$")
            phi_ky0 = np.real_if_close(fieldd.dataarray)
            erad_t = np.empty(phi_ky0.shape)
            #for tind in range(len(fieldd.timearray)):
            #    erad_t[tind, :] = -np.gradient(phi_ky0[tind, :],
            #                                   (prd.xs[1] - prd.xs[0])/prd.cm.pnt.rhostar)
            axlist[n].plot(prd.xs, av.mytrapz(erad_t, fieldd.timearray), label=r"$E_r$")
            axlist[n].set_xlabel(r'$x/a$')
            axlist[n].set_xlim(self.xmin, self.xmax)
            axlist[n].set_ylabel(r'$j_{BS}$')

    def plot_fb_param(self, axlist, prd, fieldd):
        """Plot the force balance parameter k

        Plots time averaged < u_par B >_FS for every species.
        The bootstrap current would be their sum
        :param axlist: List of axes objects (one per each species)
        :param prd: ProfileData object containing the data
        :param fieldd: FieldData object corresponding to prd
        """
        fb_param = calc_fb_param(prd, fieldd)
        for n, prs in enumerate(prd.prspec):
            axlist[n].plot(prd.xs, fb_param[n], label="GENE")
            axlist[n].plot(prd.xs, nc_predictions.fb_param_hinton_haz(prd), ":", label="H-H")
            axlist[n].plot(prd.xs, nc_predictions.fb_param_hirsch_sig(prd), "--", label="H-S")
            axlist[n].set_xlabel(r'$x/a$')
            axlist[n].set_ylabel(r'$k$')
            axlist[n].set_xlim(self.xmin, self.xmax)
            axlist[n].set_xlim(0.25, 0.85)
            axlist[n].set_ylim(0, 1.3)

    def plot_nustar(self, axlist, prd):
        """ Plot the collisionality and the inverse aspect ratio for comparison

        :param axlist: List of axes objects (one per each species)
        :param prd: ProfileData object containing the data
        """
        eps = prd.xs*prd.cm.pnt.minor_r/prd.cm.pnt.major_R
        for n, prs in enumerate(prd.prspec):
            axlist[n].plot(prd.xs, av.mytrapz(prs.nustar, prs.timearray)*eps ** 1.5,
                           label=r"$\nu_* \cdot \epsilon^{3/2}$ ")
            axlist[n].plot(prd.xs, eps ** 1.5, label=r"$\epsilon^{3/2}$")
            # B \circ Bref near the mag axis
            potwidth = nc_predictions.calc_potwidth(prd, n + 1)
            trapfr = (potwidth/prd.cm.pnt.major_R) ** 0.5
            print("Trapping fraction (potato):", trapfr)
            print("Potato width:", potwidth)
            axlist[n].hlines(trapfr ** 3, 0, potwidth)
            mspec = prd.cm.pars["mass{}".format(n + 1)]
            colltmp = []
            for tb, nu in enumerate(prs.nustar):
                colltmp.append((nu/prd.q/prd.cm.pnt.major_R*np.sqrt(
                        2*prd.prspec[0].dataarray["Ts"][tb]/mspec)*eps ** 1.5) ** (-1))
            # This is tau_ii as defined in Helander/Sigmar, note tau_i = sqrt(2) * tau_ii
            colltime = av.mytrapz(colltmp, prs.timearray)
            axlist[n].text(0.1, 0.5, "Highest collision time: {}".format(np.max(colltime)))
            axlist[n].text(0.1, 0.3, "Lowest collision time: {}".format(np.min(colltime)))
            axlist[n].set_yscale('log')
            axlist[n].set_xlabel(r'$x/a$')
            axlist[n].set_xlim(self.xmin, self.xmax)

    @classmethod
    def show(cls):
        plt.show()


class PlotProfxt(Plotting):
    """Plot contours of profiles and fluxes in x-t space.

    :param profdata: List of ProfileData objects to plot
    """

    def __init__(self, profdata):
        super().__init__()
        self.prdlist = profdata
        for prd in self.prdlist:
            if prd.isdatapresent is False:
                print("No profile data present in object ", prd.cm.fext)
        self.xmin = min([prd.xs[1] for prd in self.prdlist])
        self.xmax = max([prd.xs[-1] for prd in self.prdlist])

    def createfigs(self, poutput):
        """ Generate the figures and axes objects to plot into.

        Calls the individual plotting routines
        Currently creates one window per quantity and per species with subplots for each run
        given in the profdata list
        :param poutput: Switch to generate pdf output files
        """
        # Number of plotted variables
        numvars = 3
        # Number of plot windows
        numfigs = numvars*len(self.prdlist[0].prspec)
        # Number of subplots per window
        numsub = len(self.prdlist)
        figlist = [plt.figure() for _ in range(numfigs)]
        # One plot window for each species, but subplots for the different input runs
        axlist = [[] for _ in range(numfigs)]
        # Define a subplot grid which leaves room for the colorbar
        gs = plt.GridSpec(2, numsub, height_ratios=[12, 1])
        for ifig, fig in enumerate(figlist):
            axlist[ifig] = [fig.add_subplot(gs[0, isub]) for isub in range(numsub)]
        for iprd, prd in enumerate(self.prdlist):
            prdax = [axlist[i][iprd] for i in range(numfigs)]
            omtcm = self.plotprof_cont(prdax[0:numfigs//numvars], prd, quantity="omts")
            fluxcm = self.plotprof_cont(prdax[numfigs//numvars:2*numfigs//numvars], prd, "Qturb")
            jbscm = self.plotprof_cont(prdax[2*numfigs//numvars: 3*numfigs//numvars], prd, "jbs")
        # Generate colorbars, currently using the first set
        figlist[0].colorbar(omtcm[0], cax=figlist[0].add_subplot(gs[1, :]),
                            orientation="horizontal")
        figlist[1].colorbar(fluxcm[0], cax=figlist[1].add_subplot(gs[1, :]),
                            orientation="horizontal")
        figlist[2].colorbar(jbscm[0], cax=figlist[2].add_subplot(gs[1, :]),
                            orientation="horizontal")
        if poutput:
            figlist[0].savefig('profiles_xt{}.pdf'.format(self.prdlist[0].fext))
            figlist[1].savefig('radfluxes_xt{}.pdf'.format(self.prdlist[0].fext))
            figlist[2].savefig('jbs_xt{}.pdf'.format(self.prdlist[0].fext))

    def plotprof_cont(self, axlist, prd, quantity, vmin=0, vmax=0):
        """ Plot contour of the temperature, density or their gradient.

        :param axlist: List of axes that should have the length of the number of species
         of the second arg
        :param prd: ProfileData object which is to be plotted
        :param vmin: Minimum for the colormap
        :param vmax: Maximum for the colormap
        :param quantity: which of the fields in ProfileData to plot
        :returns: List of plotted colormaps (for colorbar)
        """
        cmlist = []
        for iprs, prs in enumerate(prd.prspec):
            axlist[iprs].set_xlim([np.min(prd.xs), np.max(prd.xs)])
            axlist[iprs].set_ylim([np.min(prs.timearray), np.max(prs.timearray)])
            colormap = self._set_colormap(quantity)
            if not (vmin == 0 and vmax == 0):
                cmlist.append(axlist[iprs].pcolormesh(prd.xs, prs.timearray,
                                                      np.array(prs.dataarray[quantity]),
                                                      shading="gouraud", cmap=colormap,
                                                      rasterized=True, vmin=vmin, vmax=vmax))
            else:
                cmlist.append(axlist[iprs].pcolormesh(prd.xs, prs.timearray,
                                                      np.array(prs.dataarray[quantity]),
                                                      shading="gouraud", cmap=colormap,
                                                      rasterized=True))
            axlist[iprs].set_title(self.titles[quantity], y=1.03)
            axlist[iprs].set_xlabel(r"$x/a$")
            axlist[iprs].set_ylabel(r"$t/(a/c_s)$")
        return cmlist

    def _set_colormap(self, quantity):
        cmapdict = {"Ts": self.cmap_unidirect, "ns": self.cmap_unidirect,
                    "omts": self.cmap_unidirect, "omns": self.cmap_unidirect,
                    "Gammanc": self.cmap_bidirect, "Gammaturb": self.cmap_bidirect,
                    "Qturb": self.cmap_unidirect, "Qnc": self.cmap_unidirect,
                    "Pinc": self.cmap_bidirect, "Piturb": self.cmap_bidirect,
                    "jbs": self.cmap_bidirect}
        return cmapdict[quantity]
