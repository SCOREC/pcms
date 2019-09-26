"""plotfsamom.py: Module to do the plotting of flux-surface averaged moments """

import numpy as np
import pydiag.utils.averages as av
from pydiag.diagplots.baseplot import plt, Plotting


class PlotFSA(Plotting):
    """ Class to plot time averaged FSA moment data

    :param fsalist: List of FSAmomdata objects to plot
    """

    def __init__(self, fsalist):
        super().__init__()
        self.fsalist = fsalist
        self.figlist = []
        self.axlist = []
        self.leglist = []
        self.labellist = ["x", "tot", r"$d/dt_{tot}$", "coll", r"$d\chi/dxy$", r"$f_1 src$",
                          r"$f_1 buff$", r"$f_0 term$", r"$v_D f_1$", r"hyp_v", r"hyp_z", r"hyp_x",
                          r"zv poisson"]

    def createfigs(self, poutput):
        """ Generate plot windows

        :param poutput: Generate the plots as pdfs as well
        """
        nfigs = len(self.fsalist[0].fsaspec)*len(self.fsalist)
        nfigsperfsa = len(self.fsalist[0].fsaspec)
        self.figlist = [plt.figure() for _ in range(nfigs)]
        # One plot window for each moment, species and quantity, one axis object in each
        self.axlist = [fig.add_subplot(111) for fig in self.figlist]
        for ifsa, fsa in enumerate(self.fsalist):
            self.plot_fsa(self.axlist[ifsa*nfigsperfsa:(ifsa+1)*nfigsperfsa], fsa)
        if poutput:
            for ifig, fig in enumerate(self.figlist):
                fig.savefig("fsamoments{}.pdf".format(ifig))

    def plot_fsa(self, axlist, fsa):
        """ Plot the results of the FSAmom diagnostic

        :param axlist: List of axes objects (one per each species)
        :param fsa: FSAmomData object containing the data
        """
        if not fsa.isdatapresent:
            print("No fsa data to plot")
            return
        Afs = (4*np.pi ** 2*self.fsalist[0].cm.pnt.major_R*self.fsalist[0].cm.pnt.minor_r /
               self.fsalist[0].cm.pnt.rhostar)
        intfsa = np.zeros((fsa.cm.pnt.nx0, fsa.fsaspec[0].fsacolumns))
        for n, fsd in enumerate(fsa.fsaspec):
            avfsa = av.mytrapz(fsd.dataarray, fsd.timearray)
            xs = fsd.dataarray[0][:, 0]
            for i in range(fsa.cm.pnt.nx0 - 1):
                intfsa[i+1, :] = intfsa[i, :] + avfsa[i+1, :]*Afs*xs[i+1]*(xs[1] - xs[0])
            intfsa = avfsa
            for iterm in range(2, 5):
                axlist[n].plot(xs, intfsa[:, iterm], label=self.labellist[iterm])
            axlist[n].legend()
            axlist[n].set_xlabel(r'$x/a$')
            axlist[n].set_title(fsd.momtitle)
