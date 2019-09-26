""" Class to do plots of the anisotropy of B """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from multiprocessing.dummy import Pool

from pydiag.diagplots.baseplot import Plotting
from pydiag.data.slices import MomFieldSlice

from pydiag.utils.comm import DiagSpace
import pydiag.utils.averages as averages
import pydiag.utils.fourier as fourier
from pydiag.utils import errors as err

class Anisotropy(object):
    """ Class for evaluating the scale-dependent anisotropy of magnetic field
    structures.
    B_L: local mean magnetic field
    b_l: local fluctuating field
    B_l: local field without filter

    :param commonlist: CommonData objects of the runs to plot
    :param rundatafiles: RunDataFiles objects of the runs to plot
    """

    def __init__(self, commonlist, rundatafiles, species=None):
        super().__init__()
        self.rundatafiles=rundatafiles
        self.cm = commonlist
        self.species = species

        #use fewer wavenumbers than ky modes for better statistics
        self.n_k=np.amin((int(2**np.rint(np.log2(self.cm.pnt.nky0/2))),32))
        print("Using {} k values for anisotropy diagnostic.".format(self.n_k))

        #calculate normalizing factor from an isotropic injection scale
        #numbers for solar wind according to Howes JGR 2008!
        k0=1e-4*np.sqrt(2) #isotropic injection scale at k*rho_i=10^-4 (rho_i defined with sqrt(2) by Howes)
        kymin=self.cm.pnt.kymin
        #assume Goldreich-Sridhar scaling of kpar~kperp^2/3 for MHD range (sda=scale-dependent anisotropy)
        #alternative: Boldyrev kpar~kperp^1/2
        model='GS'
        mhd_sda={'GS': 2./3, 'B': 1./2}
        names={'GS': 'Goldreich-Sridhar', 'B': 'Boldyrev'}
        kpar0=(kymin/k0)**(mhd_sda[model])*k0
        print("Calculating normalization factor according to {0} scaling. Result: epsilon = {1}".format(names[model],kpar0))
        #from GENE's box definitions, the normalizing epsilon factor equals kpar0
        self.epsl=kpar0

        self.nprocs=4
        self.k_array=np.zeros(self.n_k)
        self.kpar_av=np.zeros(self.n_k)
        #Cho & Lazarian: filter_factor=2
        self.filter_factor=2.

        #execute
        self.eval_kpar(species=species, rundatafiles=rundatafiles)


    def _fetch_moms(self, species, rundatafiles):
        """ Get mom and field slice object for species"""
        diagspace = DiagSpace(self.cm.spatialgrid, x_fourier=True, y_fourier=True, z_fourier=True)
        if self.cm.electromagnetic:
            moms={"apar": None}
            if self.cm.bpar:
                moms.update({"bpar": None})

        with Pool() as p:
            mom_temp = p.map(
                lambda mom: MomFieldSlice(self.cm, mom, species, diagspace, rundatafiles,
                                          modifier_func=lambda dum: dum), moms)
        for imom, mom in enumerate(moms):
            moms[mom] = mom_temp[imom]
        return moms

    def calc_individual_kpar(self,kval,a_par_fft,b_par_fft):
        k=self.k_array.tolist().index(kval)
        #print('k = {}'.format(kval))
        b_shape=np.shape(b_par_fft)
        #local mean field
        B_L=np.zeros((3,b_shape[0],b_shape[1],b_shape[2]),dtype=np.complex128)
        #local fluctuating field
        b_l=np.zeros((3,b_shape[0],b_shape[1],b_shape[2]),dtype=np.complex128)
        #aux array of ones with shape of b_l[0], used for array manipulations
        ones_b_l=np.ones(b_l[0].shape,dtype=np.float64)

        #3-d array of absolute value of k vector
        k_magn=np.sqrt(self.kx[:,None,None]**2+self.ky[None,:,None]**2+self.kz[None,None,:]**2*self.epsl**2)

        #Filters according to Cho & Lazarian 2004:
        #add to local mean field where k_magn <= kval/filter_factor
        ff=self.filter_factor
        B_L[0]=np.where(k_magn <= kval/ff, 1j*self.ky[None,:,None]*a_par_fft, 0.)
        B_L[1]=np.where(k_magn <= kval/ff, -1j*self.kx[:,None,None]*a_par_fft, 0.)
        B_L[2]=np.where(k_magn <= kval/ff, b_par_fft, 0.)

        #add to local fluctuating field where kval/2 < k_magn <= 2*kval
        b_l[0]=np.where(np.logical_and(k_magn>kval/ff,k_magn<=ff*kval), 1j*self.ky[None,:,None]*a_par_fft, 0.)
        b_l[1]=np.where(np.logical_and(k_magn>kval/ff,k_magn<=ff*kval), -1j*self.kx[:,None,None]*a_par_fft, 0.)
        b_l[2]=np.where(np.logical_and(k_magn>kval/ff,k_magn<=ff*kval), b_par_fft, 0.)

        #multiply by epsilon to transform to background units
        B_L*=self.epsl
        b_l*=self.epsl

        #compute nabla(b_l) matrix, z derivatives need extra epsilon factor to
        #compensate for different normalization
        #the grad_k step is a little awkward, but reduces the size of the code below
        grad_k=np.array([self.kx[:,None,None]*ones_b_l,
                         self.ky[None,:,None]*ones_b_l,
                         self.kz[None,None,:]*ones_b_l*self.epsl])
        grad_b_l_x=1j*grad_k*b_l[None,0,:,:,:]
        grad_b_l_y=1j*grad_k*b_l[None,1,:,:,:]
        grad_b_l_z=1j*grad_k*b_l[None,2,:,:,:]
        del ones_b_l, grad_k

        #Transform B_L, b_l, and grad(b_l) to real space for dot product
        B_L=fourier.kz_to_z(B_L,self.cm.pnt.nz0)
        b_l=fourier.kz_to_z(b_l,self.cm.pnt.nz0)
        grad_b_l_x=fourier.kz_to_z(grad_b_l_x,self.cm.pnt.nz0)
        grad_b_l_y=fourier.kz_to_z(grad_b_l_y,self.cm.pnt.nz0)
        grad_b_l_z=fourier.kz_to_z(grad_b_l_z,self.cm.pnt.nz0)

        B_L=np.fft.irfft2(B_L, axes=(-3,-2))*2*self.cm.pnt.nx0*self.cm.pnt.nky0
        b_l=np.fft.irfft2(b_l, axes=(-3,-2))*2*self.cm.pnt.nx0*self.cm.pnt.nky0
        grad_b_l_x=np.fft.irfft2(grad_b_l_x, axes=(-3,-2))*2*self.cm.pnt.nx0*self.cm.pnt.nky0
        grad_b_l_y=np.fft.irfft2(grad_b_l_y, axes=(-3,-2))*2*self.cm.pnt.nx0*self.cm.pnt.nky0
        grad_b_l_z=np.fft.irfft2(grad_b_l_z, axes=(-3,-2))*2*self.cm.pnt.nx0*self.cm.pnt.nky0

        #add background field
        B_L[2]+=1.0

        # Calculate kpar (Groselj)
        compx = np.sum(B_L*grad_b_l_x,axis=0)
        compy = np.sum(B_L*grad_b_l_y,axis=0)
        compz = np.sum(B_L*grad_b_l_z,axis=0)
        del grad_b_l_x, grad_b_l_y, grad_b_l_z

        dotproductabs = compx**2+compy**2+compz**2
        numerator=np.average(dotproductabs,axis=(-3,-2,-1))
        denom1=np.average(B_L[0]**2+B_L[1]**2+B_L[2]**2,axis=(-3,-2,-1))
        denom2=np.average(b_l[0]**2+b_l[1]**2+b_l[2]**2,axis=(-3,-2,-1))
        denominator=denom1 * denom2
        del B_L, compx, compy, compz, b_l
        return np.average(np.sqrt(numerator/denominator))

    def eval_kpar(self,species, rundatafiles):

        """ Plot contour of kpar-kperp-spectrum
        """
        #get kx,ky,kz via self.cm.spatialgrid
        self.kx=self.cm.spatialgrid.kx_fftorder
        self.ky=self.cm.spatialgrid.ky
        self.kz=self.cm.spatialgrid.kz

        #We assume an isotropic box and take ky_min and ky_max
        #as constraints for the perpendicular wavenumber range.
        #For the lowest ky modes, the background filter will be incomplete
        #and dominated by the actual background field.
        self.k_array=np.logspace(np.log10(self.ky[1]),np.log10(self.ky[-1]),self.n_k)
        ###################################
        moms = self._fetch_moms(species,rundatafiles)
        moms["apar"].check_times()
        pos = moms["apar"].calc_positions()
        self.timearray = np.take(moms["apar"].timearray, pos)

        self.kpar=np.zeros((len(self.timearray),len(self.k_array)),dtype=np.float64)
        for time in self.timearray[::max(int(len(self.timearray)/20),1)]:
            tind=self.timearray.tolist().index(time)
            print("Time: {}".format(time))

            # get b-field from field file
            b_par_fft=moms["bpar"].generate_slice_attime(time)
            a_par_fft=moms["apar"].generate_slice_attime(time)

            # Calculate Anisotropy with Fourier based filter method (Cho & Lazarian 2004)
            # Use limited number of threads, since this can easily use a whole node's memory for a large simulation
            with Pool(self.nprocs) as p:
                self.kpar[tind] = p.map(
                    lambda kval: self.calc_individual_kpar(kval, a_par_fft, b_par_fft), self.k_array)
        #compute time averaged kpar spectrum
        self.kpar_av = averages.mytrapz(self.kpar, self.timearray)*max(int(len(self.timearray)/20),1)


class PlotAnisotropy(Plotting):
    def __init__(self,commonlist,anisotropyserieslist):
        super().__init__()
        self.anisotropylist = anisotropyserieslist
        self.cm = commonlist
    def createfigs(self):
        for cm, anisotropyseries in zip(self.cm,self.anisotropylist):
            k_arr=anisotropyseries.k_array
            ff=anisotropyseries.filter_factor
            fig = plt.figure(figsize=(8.,5.))
            ax = fig.add_subplot(111)
            ax.loglog(k_arr, anisotropyseries.kpar_av,'o',label="$k_\parallel$")
            #find index closest to kperp=1
            ind=np.argmin(np.abs(k_arr-1))
            #get kpar value at kperp=1 to plot fits close to that curve
            kp_1=anisotropyseries.kpar_av[ind]
#            ax.loglog(k_arr, k_arr**(5./3)*0.1,'-',label="$k_\perp^{5/3}$")
#            ax.loglog(k_arr, k_arr**(3./3)*0.1,'-',label="$k_\perp^{3/3}$")
            ax.loglog(k_arr, k_arr**(2./3)*kp_1*1.2,'-',label="$k_\perp^{2/3}$")
#            ax.loglog(k_arr, k_arr**(1./2)*kp_1*0.8,'-',label="$k_\perp^{1/2}$")
            ax.loglog(k_arr, k_arr**(1./3)*kp_1*1.2,'--',label="$k_\perp^{1/3}$")
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            #draw shaded rectangles in wavenumber ranges where filter ranges are not completely filled with modes
            rect_l=mpatches.Rectangle((xlim[0],ylim[0]),ff*np.amin(k_arr)-xlim[0],ylim[1]-ylim[0],fc=(0.5,0.5,0.5,0.25))
            rect_r=mpatches.Rectangle((np.amax(k_arr)/ff,ylim[0]),xlim[1]-np.amax(k_arr)/ff,ylim[1]-ylim[0],fc=(0.5,0.5,0.5,0.25))
            ax.add_patch(rect_l)
            ax.add_patch(rect_r)
            ax.set_xlabel(r"$k_{\perp}\rho_i$")
            ax.set_ylabel(r"$k_{\parallel}\rho_i$")
            ax.set_title(r"Filter Method", y=1.03)
            ax.legend(loc='upper left', fontsize="medium")
            fig.tight_layout()
            fig.savefig('filter_method{}.pdf'.format(cm.fileextension))

            # name output file
#            file=open('data_kpar.txt','w')
#            file.write('#%16s %16s\n'%('kperp','kpar'))
#            for i in range(len(k_arr)):
#                file.write('%16.8e %16.8e\n'%(k_arr[i], anisotropyseries.kpar_av[i]))
#            file.close()

#    @classmethod
#    def show(cls):
#        plt.show()
