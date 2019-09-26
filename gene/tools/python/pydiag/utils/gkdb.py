""" gkdb.py: class handling data conversion and writing to the GyroKineticDataBase (https://github.com/gkdb/gkdb) """
import re
import os
import datetime
import numpy as np
import warnings
from collections import OrderedDict
import pydiag.utils.geom
import pydiag.utils.comm
import pydiag.data.omega_eigenvalue as eigenvalue
import json
import pydiag.utils.averages as averages
import matplotlib.pyplot as plt #remove later
import pydiag.diagplots.plot_ball as ballooning
import pydiag.data.nrgdata as nrg

class GKDB(object):
    """ GKDB class:

    Converts the results of a GENE run into a GKDB-style json file

    Currently restricted to rotationless linear fluxtube simulations with Miller equilibria"
    """

    def __init__(self,common,rundatafiles):
        # dictionary for gkdb parameters: {parameter: value}
        self.cm = common
        self.geom = pydiag.utils.geom.Geometry(self.cm)
        self.fields = ["phi"]
        if common.electromagnetic:
            self.fields += ["apar"]
            if common.bpar:
                self.fields += ["bpar"]

        self.rundatafiles = rundatafiles
        self.gkdict = OrderedDict()
        self.set_gkdb_shape()
        self.set_gene2gkdb_conversion()
        self.fill_gkdb_from_par()

    def set_gkdb_shape(self):
        """ Determines shape parameters for GKDB """

        print ("Determining flux surface shape ...")

        nz0 = self.cm.pnt.nz0
        self.zval = 2.0*np.pi*(-0.5+list(range(nz0))/np.real(nz0))

        #determine flux surface center
        self.R0 = 0.5*(max(self.geom.R)+min(self.geom.R))
        self.Z0 = 0.5*(max(self.geom.Z)+min(self.geom.Z))
        print ("(R0,Z0) = (", self.R0, ",", self.Z0,")")
#        self.R0 = np.sum(self.geom.R*self.geom.jacobian)/np.sum(self.geom.jacobian)
#        self.Z0 = np.sum(self.geom.Z*self.geom.jacobian)/np.sum(self.geom.jacobian)

        #minor radius of flux surface
        self.r = 0.5*(max(self.geom.R)-min(self.geom.R))

        #poloidal theta (COCOS=11)
        self.theta = np.arctan(-(self.geom.Z-self.Z0)/(self.geom.R-self.R0))
        for iz in range(int(nz0/2),nz0-1):
            thetadiff = (self.theta[iz+1]-self.theta[iz])
            if thetadiff>0.5*np.pi:
                self.theta[iz+1]-=np.pi
            elif thetadiff<-0.5*np.pi:
                self.theta[iz+1]+=np.pi
            thetadiff = (self.theta[nz0-1-iz-1]-self.theta[nz0-1-iz])
            if thetadiff>0.5*np.pi:
                self.theta[nz0-1-iz-1]-=np.pi
            elif thetadiff<-0.5*np.pi:
                self.theta[nz0-1-iz-1]+=np.pi

        #sorted theta (starting at th=0)
        theta_sort = (self.theta+2.0*np.pi)%(2.0*np.pi)
        sort_ind = np.argsort(theta_sort)
        theta_sort = theta_sort[sort_ind]

        #Flux surface distance (original grid and normalized)
        aN = np.sqrt((self.geom.R[sort_ind]-self.R0)**2+(self.geom.Z[sort_ind]-self.Z0)**2)/self.R0
        aN_dr = ((self.geom.R[sort_ind]-self.R0)*self.geom.dxdR[sort_ind]+
                 (self.geom.Z[sort_ind]-self.Z0)*self.geom.dxdZ[sort_ind])/\
            (aN*self.geom.gxx[sort_ind])/self.R0

        #equidistant theta (COCOS=11) with somewhat higher resolution for fourier transforms etc
        nth = 4*nz0
        theta_grid = (np.arange(0.0,nth,1.0))*2.0*np.pi/nth

        #interpolate to this equidist. grid with proper order (interpol requires increasing arr)
        aN_th = np.interp(theta_grid,theta_sort,aN,period=2.0*np.pi)
        aN_dr_th = np.interp(theta_grid,theta_sort,aN_dr,period=2.0*np.pi)

        #now fourier transform
        fft_aN = np.fft.fft(aN_th)
        fft_aN_dr = np.fft.fft(aN_dr_th)

        #keep up to nth_kept modes (typically around 8)
        nth_kept = 8 #nz0
        self.cN = 2.0*np.real(fft_aN)[0:nth_kept]/nth
        self.cN[0] *= 0.5
        self.sN = -2.0*np.imag(fft_aN)[0:nth_kept]/nth
        self.sN[0] *= 0.5
        self.cN_dr = 2.0*np.real(fft_aN_dr)[0:nth_kept]/nth
        self.cN_dr[0] *= 0.5
        self.sN_dr = -2.0*np.imag(fft_aN_dr)[0:nth_kept]/nth
        self.sN_dr[0] *= 0.5

        #check parametrization
        nind = np.arange(0.0,nth_kept,1.0)
        a_check = []
        a_dr_check = []
        for iz in range(0,nth):
            a_check += [np.sum(self.cN*np.cos(nind*theta_grid[iz])+
                               self.sN*np.sin(nind*theta_grid[iz]))]
            a_dr_check += [np.sum(self.cN_dr*np.cos(nind*theta_grid[iz])+
                                  self.sN_dr*np.sin(nind*theta_grid[iz]))]

        a_dr_check = np.array(a_dr_check)
        err = np.sum(np.abs(a_check-aN_th)/aN_th)/nth
        print ("Relative error of flux surface parametrization: {0:12.6E}".format(err))

        #some plot diagnostics
        show_aN = False
        if show_aN:
            plt.plot(theta_sort,aN,label='aN_orig')
            plt.plot(theta_grid,aN_th,label='aN_th,intermediate')
            plt.plot(theta_grid,a_check,label='aN_final')
            plt.plot(theta_sort,aN_dr,label='aN_dr_orig')
            plt.plot(theta_grid,a_dr_check,label='aN_dr_final')
            plt.legend()
            plt.show()

        show_flux_surface = True
        if show_flux_surface:
            dr = 0.04  #radial extension for testing

            plt.plot(self.geom.R[sort_ind],self.geom.Z[sort_ind],'o',label="original")
            plt.plot((self.geom.R[sort_ind]+dr*self.geom.dxdR[sort_ind]/self.geom.gxx[sort_ind]),
                     (self.geom.Z[sort_ind]+dr*self.geom.dxdZ[sort_ind]/self.geom.gxx[sort_ind]),
                     'x',label="original+dr")
            theta_grid = (np.arange(0.0,nth,1.0))*2.0*np.pi/(nth-1)
            # R_n_th = self.R0+aN_th*np.cos(theta_grid)
            # Z_n_th = self.Z0-aN_th*np.sin(theta_grid)
            # plt.plot(R_n_th,Z_n_th,label="intermediate")
            R_check = self.R0+a_check*np.cos(theta_grid)*self.R0
            Z_check = self.Z0-a_check*np.sin(theta_grid)*self.R0
            plt.plot(R_check,Z_check,label="gkdb")
            R_check = self.R0+(a_check+dr/self.R0*a_dr_check)*np.cos(theta_grid)*self.R0
            Z_check = self.Z0-(a_check+dr/self.R0*a_dr_check)*np.sin(theta_grid)*self.R0
            plt.plot(R_check,Z_check,label="gkdb+dr")
            plt.xlabel('R/m')
            plt.ylabel('Z/m')
            plt.axis('equal')
            plt.legend()
            plt.show()

    def set_gene2gkdb_conversion(self):
        """ set conversion factors for various quantities """

        print ("Setting conversion factors ...")

        #determine index of electron species
        ispec_electrons = -1111
        for spec in self.cm.specnames:
            ispec = self.cm.specnames.index(spec) + 1
            if self.cm.pars["charge{}".format(ispec)]==-1:
                ispec_electrons = ispec

        self.gene2gkdb = OrderedDict()
        self.gene2gkdb["Lref"] = self.cm.pnt.Lref/self.R0 #Lref_GENE/Lref_GKDB
        self.gene2gkdb["Bref"] = 1.0  #FILL
        self.gene2gkdb["dxdr"] = 1.0  #FILL
        self.gene2gkdb["Tref"] = 1.0/self.cm.pars["temp{}".format(ispec_electrons)]
        self.gene2gkdb["mref"] = self.cm.pnt.mref/1.999007501778479 #mref_GENE/mref_GKDB (in units of proton mass)
        self.gene2gkdb["nref"] = 1.0/self.cm.pars["dens{}".format(ispec_electrons)]
        self.gene2gkdb["qref"] = 1.0

        print ("Lref_GENE/Lref_GKDB = ", self.gene2gkdb["Lref"])

        #derived quantity conversion
        self.gene2gkdb["vth"]  = np.sqrt(0.5*self.gene2gkdb["Tref"]/self.gene2gkdb["mref"]) #cref_GENE/vth_GKDB
        self.gene2gkdb["rho"]  = self.gene2gkdb["mref"]*self.gene2gkdb["vth"]/(self.gene2gkdb["qref"]*self.gene2gkdb["Bref"])
        self.gene2gkdb["rho_Lref"] = self.gene2gkdb["rho"]/self.gene2gkdb["Lref"]

        print ("rho_Lref_GENE/rho_Lref_GKDB = ", self.gene2gkdb["rho_Lref"])

        #field conversion
        self.gene2gkdb["phi"] = self.gene2gkdb["Tref"]/self.gene2gkdb["qref"]*self.gene2gkdb["rho_Lref"]
        self.gene2gkdb["apar"] = self.gene2gkdb["Lref"]*self.gene2gkdb["Bref"]*self.gene2gkdb["rho_Lref"]**2
        self.gene2gkdb["bpar"] = self.gene2gkdb["Bref"]*self.gene2gkdb["rho_Lref"]

        #moments conversions (may be species dependent)
        self.gene2gkdb["dens"] = []
        self.gene2gkdb["upar"] = []
        self.gene2gkdb["tpar"] = []
        self.gene2gkdb["tperp"] = []
        for spec in self.cm.specnames:
            ispec = self.cm.specnames.index(spec) + 1
            self.gene2gkdb["dens"] += [self.gene2gkdb["nref"]*self.gene2gkdb["rho_Lref"]*\
                                       self.cm.pars["dens{}".format(ispec)]]
            self.gene2gkdb["upar"] += [self.gene2gkdb["vth"]*self.gene2gkdb["rho_Lref"]]
            self.gene2gkdb["tpar"] += [self.gene2gkdb["Tref"]*self.gene2gkdb["rho_Lref"]*\
                                       self.cm.pars["temp{}".format(ispec)]]
            self.gene2gkdb["tperp"] += [self.gene2gkdb["Tref"]*self.gene2gkdb["rho_Lref"]*\
                                        self.cm.pars["temp{}".format(ispec)]]

        #fluxes
        self.gene2gkdb["pflux"] = self.gene2gkdb["nref"]*self.gene2gkdb["vth"]*self.gene2gkdb["rho_Lref"]**2
        self.gene2gkdb["mflux"] = self.gene2gkdb["nref"]*self.gene2gkdb["mref"]*self.gene2gkdb["Lref"]*\
                                  (self.gene2gkdb["vth"]*self.gene2gkdb["rho_Lref"])**2   #LREF?
        self.gene2gkdb["qflux"] = self.gene2gkdb["Tref"]*self.gene2gkdb["pflux"]

    def fill_gkdb_from_par(self):
        """ Fill gkdb dict """
        pars = self.cm.pars
        geom = self.geom
        rat = self.gene2gkdb

        self.gkdict.clear()
        self.gkdict["point"] = {"creator" : os.environ['USER'],
                                "date" : str(datetime.date.today()),
                                "comment" : "json entry for testing only - note that "+\
                                "mode structure moments and flux-surface averaged fluxes "+\
                                "are given in particle space (i.e. incl. the full pullback operator)\n "+\
                                "furthermore, the index for the mode structure moments is adjusted to "+\
                                "[N_pol,grid x Nspecies x Nmodes] which seems to be the more natural order\n "+\
                                "finally, EM flux contributions are for now not separated and all attributed "+\
                                "to the Apar term"
                            }
        self.gkdict["code"] = {"name": "GENE",
                               "version" : pars['RELEASE']+' - '+(pars['GIT_MASTER'])[0:9],
#                               "parameters" : pars,
                               "include_centrifugal_effects" : pars['with_centrifugal'],
                               "include_a_parallel" : self.cm.electromagnetic,
                               "include_b_field_parallel" : self.cm.bpar,
                               "collision_pitch_only" : pars['collision_op']=="'pitch_angle'",
                               "collision_ei_only" : False,
                               "collision_momentum_conservation" : (pars['coll_cons_model']!='none'
                                                                    if pars['coll'] > 0.0 else False),
                               "collision_energy_conservation" :  (pars['coll_cons_model']!='none'
                                                                   if pars['coll'] > 0.0 else False),
                               "collision_finite_larmor_radius" : pars['collision_op']=='sugama', #check this
                               "initial_value_run": pars['comp_type']=="'IV'"}

        s_Ip = -pars['sign_Ip_CW']
        s_Bt = -pars['sign_Bt_CW']

        self.gkdict["flux_surface"] = {
            "r_minor" : pars['x0'], #CHECK
            # Derived from Shape
            "elongation" : (pars['kappa'] if pars['magn_geometry']=="'miller'" else None),
            "triangularity" : (pars['delta'] if pars['magn_geometry']=="'miller'" else None),
            "squareness" : (pars['zeta'] if pars['magn_geometry']=="'miller'" else None),
            # Non-derived
            "q" : s_Ip*s_Bt*np.abs(pars['q0']),
            "magnetic_shear" : pars['shat'], #add factor r/x dx/dr
            "pressure_gradient" : pars['dpdx_pm']*rat["Bref"]**2*rat["dxdr"],
            "ip_sign" : s_Ip,
            "b_field_tor_sign" : s_Bt,
            # Original shape
            "c" : self.cN.tolist(), #NORMALIZATION?
            "s" : self.sN.tolist(), #NORMALIZATION?
            "dc_dr_minor" : self.cN_dr.tolist(),
            "ds_dr_minor" : self.sN_dr.tolist()}

        coulomb_log = 24.-np.log(np.sqrt(self.cm.pnt.nref*1E13)/
                                self.cm.pnt.Tref*0.001)

        self.gkdict["species_global"] = {
            "beta" : pars['beta']*rat["nref"]*rat["Tref"]*rat["Bref"]**2,
            "collisionality" : pars['coll']*8.0*np.sqrt(2.0)*rat["Tref"]**2/\
            (rat["qref"]**4*rat["Lref"]*rat["nref"])/coulomb_log,
            "collision_enhancement_factor" : 1.0,
            "toroidal_velocity" : pars["Omega0_tor"]*rat["vth"]/rat["Lref"],
            "debye_length" : np.sqrt(pars['debye2']), #ADD FACTORS
            # Derived from Species
            "zeff" : 1.0 #pars['Zeff'] #self-consistent computation?
        }

        self.gkdict["species"] = []
        for ispec in range(1,pars['n_spec']+1):
            self.gkdict["species"] += [{
                "name" : self.cm.pars["name{}".format(ispec)],
                "charge" : self.cm.pars["charge{}".format(ispec)],
                "mass" : self.cm.pars["mass{}".format(ispec)]*rat["mref"],
                "density" : self.cm.pars["dens{}".format(ispec)]*rat["nref"], #zero if passive species?
                "temperature" : self.cm.pars["temp{}".format(ispec)]*rat["Tref"],
                "density_log_gradient" : self.cm.pars["omn{}".format(ispec)]/rat["Lref"]*rat["dxdr"],
                "temperature_log_gradient" : self.cm.pars["omt{}".format(ispec)]/rat["Lref"]*rat["dxdr"],
                "toroidal_velocity_gradient" : pars["ExBrate"]*pars["q0"]/pars["x0"]*\
                rat["vth"]/rat["Lref"]*rat["dxdr"]
            }]

        Cyq0_x0 = geom.Cy * self.cm.pnt.q0 / self.cm.pnt.x0
        if self.cm.pnt.adapt_lx:
            absnexc = 1
        else:
            absnexc = int(np.round(self.cm.pnt.lx*self.cm.pnt.n_pol*
                    np.abs(self.cm.pnt.shat)*self.cm.pnt.kymin*
                    np.abs(Cyq0_x0)))
        nconn = int(int(int(self.cm.pnt.nx0-1)/2)/(absnexc*self.cm.pnt.ky0_ind))*2+1

        #transformation of wavevector components on a z grid
        #following D. Told, PhD thesis, Sec. A.3, p. 167
        k1facx = np.sqrt(geom.gxx)
        k1facy = geom.gxy/np.sqrt(geom.gxx)
        k2fac = np.sqrt(geom.gyy-geom.gxy**2/geom.gxx)

        #now interpolate to theta = 0
        #interpolate to equidist. grid with proper order (interpol requires increasing arr)
        sort_ind = np.argsort(self.theta)
        k1facx_th0 = np.interp(0.0,self.theta[sort_ind],k1facx[sort_ind])/rat["rho"]
        k1facy_th0 = np.interp(0.0,self.theta[sort_ind],k1facy[sort_ind])/rat["rho"]
        k2fac_th0 = np.interp(0.0,self.theta[sort_ind],k2fac[sort_ind])/rat["rho"]

        ky = self.cm.pnt.kymin*self.cm.pnt.ky0_ind
        k1val = (self.cm.pnt.kx_center*k1facx_th0+ky*k1facy_th0)
        k2val = ky*k2fac_th0

        print ("k1facx_th0, k1facy_th0, k2fac_th0: ", k1facx_th0, k1facy_th0, k2fac_th0)


        fieldmodestruct = ballooning.ModeStructure(self.cm, self.rundatafiles)
        theta_full = []
        for iconn in range(-int(nconn/2),int(nconn/2)+1):
            theta_full += list(self.theta-np.sign(self.theta[0])*iconn*2.0*np.pi)
        fullsort_ind = np.argsort(theta_full)
        thetasort = np.array(theta_full)[fullsort_ind]

        #determine normalization amplitude and phase for parallel mode structures
        Af_gkdb = 0.0
        for field in self.fields:
            Af_gkdb += np.trapz(np.abs(fieldmodestruct.ballamps[field][fullsort_ind]*rat[field])**2,thetasort)
        Af_gkdb = np.sqrt(Af_gkdb)/(2.0*np.pi) #np.sqrt((2.0*np.pi))
        alpha_gkdb = -np.interp(0.0,thetasort,np.angle(fieldmodestruct.ballamps["phi"][fullsort_ind]))
        print ("Af_gkdb = ", Af_gkdb, ", alpha_gkdb = ", alpha_gkdb)

        var_check = "phi"
        amp_check = np.sqrt(np.trapz(np.abs(fieldmodestruct.ballamps[var_check][fullsort_ind]*
                                            rat[var_check]*np.exp(1j*alpha_gkdb))**2,thetasort))/(2.0*np.pi)/Af_gkdb
        print ("amp_check = ", amp_check)

        newfieldnames = {"phi" : "phi_potential_perturbed", "apar" : "a_parallel_perturbed",
                    "bpar" : "b_field_parallel_perturbed"}

        eigendat = eigenvalue.Eigenvaluedata(self.cm)
        eigenvalues = [{
            "growth_rate" : eigendat.growth_rate*rat["vth"]/rat["Lref"],
            "frequency" : eigendat.frequency*rat["vth"]/rat["Lref"],
            "growth_rate_tolerance": self.cm.pnt.omega_prec*rat["vth"]/rat["Lref"]
            }]

        eigenvector = {"poloidal_angle" : thetasort.tolist()}
        for var in fieldmodestruct.ballamps.keys():
            print ("writing mode structure of ", var)
            eigenvector["r_"+newfieldnames[var]]=np.real(fieldmodestruct.ballamps[var][fullsort_ind]*
                                                         rat[var]*np.exp(1j*alpha_gkdb)/Af_gkdb).tolist()
            eigenvector["i_"+newfieldnames[var]]=np.imag(fieldmodestruct.ballamps[var][fullsort_ind]*
                                                         rat[var]*np.exp(1j*alpha_gkdb)/Af_gkdb).tolist()

        newmomnames = { "dens" : "density_moment",
                        "upar" : "parallel_velocity_moment",
                        "tpar" : "parallel_temperature_moment",
                        "tperp" : "perpendicular_temperature_moment"}

        for var in newmomnames.keys():
            eigenvector["r_"+newmomnames[var]] = []
            eigenvector["i_"+newmomnames[var]] = []

        for spec in self.cm.specnames:
            ispec = self.cm.specnames.index(spec) #w/o +1
            mommodestruct = ballooning.ModeStructure(self.cm, self.rundatafiles, moms=list(newmomnames.keys()), species=spec)
            for var in newmomnames.keys():
                eigenvector["r_"+newmomnames[var]] += [np.real(mommodestruct.ballamps[var][fullsort_ind]*
                                                               rat[var][ispec]*np.exp(1j*alpha_gkdb)/Af_gkdb).tolist()]
                eigenvector["i_"+newmomnames[var]] += [np.imag(mommodestruct.ballamps[var][fullsort_ind]*
                                                               rat[var][ispec]*np.exp(1j*alpha_gkdb)/Af_gkdb).tolist()]


#        eigenvalues += [{"eigenvector" : eigenvector}]
        eigenvalues[0]["eigenvector"] = eigenvector
        self.gkdict["wavevectors"] = [{
            "radial_wavevector" : k1val,
            "binormal_wavevector" : k2val,
            "poloidal_turns" : nconn,
            "eigenvalues" : eigenvalues
        }]


        nrgdat = self.rundatafiles.get_fileobject("nrg")
        nrgdat.generate_timeseries()

        self.gkdict["particle_fluxes"] = {
            "phi_potential": [(np.swapaxes(nrgdat.dataarray[:,:,4],0,1)*rat["pflux"]*(1.0/Af_gkdb)**2).tolist()],
            "a_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,5],0,1)*rat["pflux"]*(rat["rho_Lref"]/Af_gkdb)**2).tolist()],       #FIX THIS
            "b_field_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,5],0,1)*rat["pflux"]*(1.0/Af_gkdb)**2*0.0).tolist()]  #FIX THIS -- SET TO ZERO
        }
        self.gkdict["momentum_fluxes_rotating"] = {
            "phi_potential": [(np.swapaxes(nrgdat.dataarray[:,:,8],0,1)*rat["mflux"]*(1.0/Af_gkdb)**2).tolist()],
            "a_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,9],0,1)*rat["mflux"]*(1.0/Af_gkdb)**2).tolist()],       #FIX THIS
            "b_field_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,9],0,1)*rat["mflux"]*(1.0/Af_gkdb)**2*0.0).tolist()]  #FIX THIS -- SET TO ZERO
        }
        self.gkdict["heat_fluxes_rotating"] = {
            "phi_potential": [(np.swapaxes(nrgdat.dataarray[:,:,6],0,1)*rat["qflux"]*(1.0/Af_gkdb)**2).tolist()],
            "a_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,7],0,1)*rat["qflux"]*(1.0/Af_gkdb)**2).tolist()],       #FIX THIS
            "b_field_parallel": [(np.swapaxes(nrgdat.dataarray[:,:,7],0,1)*rat["qflux"]*(1.0/Af_gkdb)**2*0.0).tolist()]  #FIX THIS -- SET TO ZERO
        }


    def Write_json(self, path):
        """ Take the dict and write a gkdb json parameters file """
#        print(json.dumps(self.gkdict, indent=4))
        with open('data.json', 'w') as fp:
            json.dump(self.gkdict, fp, indent=4,separators=(',', ': '))
            #json.dump(self.gkdict, fp, sort_keys=True,indent=4,separators=(',', ': '))
