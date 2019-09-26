""" Classes for plotting field based mode structures """

import matplotlib.pyplot as plt
import numpy as np
import warnings
from multiprocessing.dummy import Pool

from pydiag.diagplots.baseplot import Plotting
from pydiag.data.slices import MomFieldSlice
from pydiag.utils.comm import DiagSpace

import pydiag.utils.averages as averages
import pydiag.utils.geom


class ModeStructure(object):
    """ Ballooning Mode Structure for selected ky of field data """

    def __init__(self, common, rundatafiles, species=None, moms=None, kyind=None, normalize=False):
        if not common.x_local or not common.y_local:
            raise NotImplementedError("Not yet implemented for nonlocal simulations")
        self.cm = common
        self.species = species
        self.timearray = None
        self.normalize = normalize

        if moms is None:
            self.ballamps = {"phi": []}
            if common.electromagnetic:
                self.ballamps.update({"apar": []})
                if common.bpar:
                    self.ballamps.update({"bpar": []})
        else:
            self.ballamps = {}
            for mom in moms:
                self.ballamps.update({mom : []})

        if kyind is None:
            self.kyind=(1 if self.cm.pnt.nonlinear else 0) #default: 1st finite ky
        else:
            self.kyind=kyind[0]

        self.calc_modestructure(species=species,rundatafiles=rundatafiles)

    def _fetch_moms(self, rundatafiles, species):
        """ Get mom and field slice object for species"""
        diagspace = DiagSpace(self.cm.spatialgrid, x_fourier=True, y_fourier=True)

        momamps = {}
        for mom in self.ballamps.keys():
            momamps.update({mom : None})

        with Pool() as p:
            momamp_temp = p.map(
                lambda mom: MomFieldSlice(self.cm, mom, species, diagspace, rundatafiles,
                                          modifier_func=lambda dum: dum), momamps)
        for imom, mom in enumerate(momamps):
            momamps[mom] = momamp_temp[imom]
        return momamps

    def calc_modestructure(self,species,rundatafiles):
        """ Calculate the mode structures"""
        momamps = self._fetch_moms(rundatafiles, species=species)
        geom = pydiag.utils.geom.Geometry(self.cm)
        timerefvar = list(self.ballamps.keys())[0]  #TODO: NEED TO SET TO LOWEST TIME RES. VAR/FILE
        momamps[timerefvar].check_times()
        pos = momamps[timerefvar].calc_positions()
        self.timearray = np.take(momamps[timerefvar].timearray, pos)

        if self.cm.pnt.magn_geometry == "'s_alpha'" or self.cm.pnt.magn_geometry == "'circular'":
            Cyq0_x0 = 1.0
        else:
            Cyq0_x0 = geom.Cy * self.cm.pnt.q0 / self.cm.pnt.x0
        sign_shear = (-1 if self.cm.pnt.shat < 0 else 1)
        nexc_sign = sign_shear*self.cm.pnt.sign_Ip_CW*self.cm.pnt.sign_Bt_CW

        if self.cm.pnt.adapt_lx:
            nexc = 1
        else:
            nexc = int(np.round(self.cm.pnt.lx*self.cm.pnt.n_pol*
                    np.absolute(self.cm.pnt.shat)*self.cm.pnt.kymin*
                    np.absolute(Cyq0_x0)))
        nexc*=nexc_sign
        nexcj = nexc*(self.cm.pnt.ky0_ind+self.kyind)
        nconn = int(int(int(self.cm.pnt.nx0-1)/2)/np.abs(nexcj))*2+1
        phasefac = (-1)**nexcj*np.exp(2.0*np.pi*1j*self.cm.pnt.n0_global*
                                         self.cm.pnt.q0*(self.cm.pnt.ky0_ind+self.kyind))

        self.chi_ext = []
        for ncind in range(-int(nconn/2),int(nconn/2)+1):
            self.chi_ext = np.concatenate([self.chi_ext,
                                           (self.cm.spatialgrid.z+ncind*2.0*np.pi)])
        for time in self.timearray:
            def mom_t(momname):
                return momamps[momname].generate_slice_attime(time)
            for mom in momamps:
                if self.normalize:
                    normval = mom_t(mom)[0,self.kyind,int(self.cm.pnt.nz0/2)]
                else:
                    normval = 0.
                balloon = []
                for ncind in range(-int(nconn/2),int(nconn/2)+1):
                    balloon = np.concatenate([balloon,np.conj(phasefac)**ncind*
                                              mom_t(mom)[ncind*nexcj,self.kyind]])
                if normval == 0.:
                    normphase=0.
                    normval=1.
                else:
                    normphase = np.angle(normval)
                    normval=np.absolute(normval)
                self.ballamps[mom].append(balloon*np.exp(-1j*normphase)/normval)

        for mom in self.ballamps:
            self.ballamps[mom] = averages.mytrapz(self.ballamps[mom], self.timearray)


class PlotModeStructure(Plotting):

    def __init__(self, commonlist, ballserieslist):
        super().__init__()
        self.balllist = ballserieslist
        self.cm = commonlist
        self.titles.update({"phi": r"$\phi$", "apar": r"$A_\parallel$", "bpar": r"$B_\parallel$"})

    def createfigs(self):
        for ballseries in self.balllist:
            moms = ballseries.ballamps.keys()
            fig = plt.figure("Ballooning Mode Structure")

            chi_pi = ballseries.chi_ext/np.pi
            kyval = ballseries.cm.spatialgrid.ky[ballseries.kyind]

            for imom,mom in enumerate(moms):
                ax_lin = fig.add_subplot(1,len(moms),imom+1)
                ballarr = ballseries.ballamps[mom]
                ax_lin.plot(chi_pi, np.absolute(ballarr),
                            label="|{}|".format(self.titles[mom]),
                            color='black')
                ax_lin.plot(chi_pi, np.real(ballarr),
                            label="Re({})".format(self.titles[mom]),
                            color='blue',lw=1)
                ax_lin.plot(chi_pi, np.imag(ballarr),
                            label="Im({})".format(self.titles[mom]),
                            color='red',lw=1)
                ax_lin.set_xlabel(r"$\chi/\pi$")
                ax_lin.set_title(r"{0:s}$(k_y\rho={1:6.3f})$".
                                 format(self.titles[mom],kyval))
                ax_lin.axhline(y=0,color='black',lw=1,ls='--')
                ax_lin.legend()

            fig.tight_layout()
