""" Module for plotting of x-y contours"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import warnings
import multiprocessing

import pydiag.diagplots.baseplot
from pydiag.data.slices import MomFieldSlice
from pydiag.utils.comm import DiagSpace

import pydiag.utils.averages as averages
import pydiag.utils.geom


class MomFieldContour(object):
    """ Create a representation of a x-y map for a mom or field quantity"""

    def __init__(self, common, species, rundatafiles, zind=None, zavg=False, quantities=("phi",)):
        if quantities is None:
            raise RuntimeError("No quantities given for contours")
        if not zind:
            zind = int(common.pnt.nz0/2)
        zrange = (zind, zind + 1)
        diagspace = DiagSpace(common.spatialgrid, False, False, False, (None, None), (None, None),
                              zrange, False, False, zavg)
        self.species = species
        self.cm = common
        self.contourdata = {}
        self.fetch_data(diagspace, quantities, rundatafiles)

    def fetch_data(self, diagspace, quantities, rundatafiles):
        if set(quantities).intersection({"phi", "apar", "bpar"}) and set(quantities).intersection(
                {"dens", "tpar", "tperp", "qpar", "qperp", "upar"}):
            bothmomandfield = True
        else:
            bothmomandfield = False
        for quantity in quantities:
            qslice = MomFieldSlice(self.cm, quantity, species=self.species, diagspace=diagspace,
                                   rundatafiles=rundatafiles, modifier_func=lambda dum: dum)

            if bothmomandfield and quantity in ["phi", "apar", "bpar"]:
                sparse_field = int(self.cm.pnt.istep_mom/self.cm.pnt.istep_field)
                qslice.generate_timeseries(sparsefactor=sparse_field)
            else:
                qslice.generate_timeseries()
            self.contourdata.update({quantity: qslice})


class PlotXYContour(pydiag.diagplots.baseplot.Plotting):
    """ Class to handle contour plots of mom or field data in x-y real space"""

    def __init__(self, contourlist):
        super().__init__()
        self.contourlist = contourlist
        self.anims = []  # The animation objects need to be stored for the animations to run
        self.figlist = []
        self.axlist = []
        self.axcblist = []

    def createfigs(self, animate=False, save=False):
        for cont in self.contourlist:
            for ic, contdata in enumerate(cont.contourdata):
                if not animate:
                    cont.contourdata[contdata].generate_timeaverage()
                    absval = int(
                            self.round_to_n(np.max(np.abs(cont.contourdata[contdata].timeaverage)),
                                            2))
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    cm1 = ax.contourf(cont.cm.spatialgrid.x, cont.cm.spatialgrid.y,
                                      np.squeeze(cont.contourdata[contdata].timeaverage).T, 100,
                                      cmap=self.cmap_bidirect, vmax=absval, vmin=-absval,
                                      zorder=-20)
                    ax.set_rasterization_zorder(z=-10)
                    ax.set_xlabel(r"$x/\rho_{ref}$")
                    ax.set_ylabel(r"$y/\rho_{ref}$")
                    ax.set_title(self.titles[contdata])
                    fig.colorbar(cm1)
                    fig.tight_layout()
                    if save:
                        fig.savefig("contour_{}{}_{}.pdf".format(contdata, cont.species,
                                                                 cont.cm.fileextension))
                else:
                    absval = int(
                        self.round_to_n(np.max(np.abs(cont.contourdata[contdata].dataarray)), 2))
                    self.figlist.append(plt.figure())
                    self.axlist.append(self.figlist[-1].add_subplot(111))

                    cma = self.axlist[-1].contourf(cont.cm.spatialgrid.x, cont.cm.spatialgrid.y,
                                                   np.squeeze(cont.contourdata[contdata].dataarray[
                                                                  0]).T, 100,
                                                   cmap=self.cmap_bidirect, vmax=absval,
                                                   vmin=-absval, zorder=-20)
                    self.axcblist.append(self.figlist[-1].colorbar(cma).ax)

                    def on_click(event):
                        if event.dblclick:
                            for ani in self.anims:
                                if ani.running:
                                    ani.event_source.stop()
                                else:
                                    ani.event_source.start()
                                ani.running ^= True

                    def animate(t, axs=None, fig=None, axcb=None, contdata=None):
                        axs.clear()
                        axcb.clear()
                        axs.set_xlabel(r"$x/\rho_{ref}$")
                        axs.set_ylabel(r"$y/\rho_{ref}$")
                        axs.set_rasterization_zorder(z=-10)
                        axs.set_title(r"{}, t={:.3f}".format(self.titles[contdata],
                                                         cont.contourdata[contdata].timearray[t]))
                        cma = axs.contourf(cont.cm.spatialgrid.x, cont.cm.spatialgrid.y,
                                           np.squeeze(cont.contourdata[contdata].dataarray[t]).T,
                                           100, cmap=self.cmap_bidirect, vmax=absval, vmin=-absval,
                                           zorder=-20)
                        cb = fig.colorbar(cma, cax=axcb)
                        return cma, cb

                    self.anims.append(animation.FuncAnimation(self.figlist[-1], animate, range(
                            len(cont.contourdata[contdata].timearray)), fargs=(
                        self.axlist[ic], self.figlist[ic], self.axcblist[ic], contdata),
                                                              blit=False))
                    for anim in self.anims:
                        anim.running = True
                    for fig in self.figlist:
                        fig.canvas.mpl_connect('button_press_event', on_click)
                    if save:
                        # Keep the plotting from freezing by doing the saving in a subprocess
                        def save_anim():
                            print("Saving animation in the background, this can take a long while")
                            self.anims[-1].save(
                                    "animatedcontour_{}_{}{}.mp4".format(contdata, cont.species,
                                                                         cont.cm.fileextension),
                                    dpi=300)
                            print("Finished saving animation")
                        p = multiprocessing.Process(target=save_anim)
                        p.daemon = True
                        p.start()
