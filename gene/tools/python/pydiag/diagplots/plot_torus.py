import numpy as np
import warnings
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.animation as animation

import pydiag.utils.comm
import pydiag.utils.geom
import pydiag.utils.averages
import pydiag.data.slices
import pydiag.diagplots.baseplot

import multiprocessing


class TorusContourCut(pydiag.data.slices.MomFieldSlice):
    """ Generate a phi (toroidal angle) cut of mom or field data and a triangulation """

    def __init__(self, common, species, rundatafiles, torangle=0, z_res=400, quantity="phi"):
        """ Constructor for torus diagnostic specific slice

        :param common: A CommonData object for the run to diagnose
        :param species: Name of the species to diagnose
        :param rundatafiles: RunDataFiles object, matching common
        :param torangle: The toroidal angle of the cut, phi=[0, 2*pi]
        :param z_res: The number of points that z should be interpolated to
        :param quantity: Name of the quantity to plot (from field or mom, refer to MomFieldSlice
        constructor
        """
        diagspace = pydiag.utils.comm.DiagSpace(common.spatialgrid, False, False, False,
                                                (None, None), (None, None), (None, None), False,
                                                False, False)
        super().__init__(common, quantity, species, diagspace, rundatafiles, lambda dum: dum)
        if common.pnt.magn_geometry == "miller":
            warnings.warn(
                    "TorusContourCut class lacks a correction factor for Miller geomtery")  # See
            #  Cyq0_r0 in IDL diag torus.pro

        self.torangle = torangle%(2*np.pi)
        self.species = species
        self.nz_extended = z_res
        self.z_full = np.append(self.cm.spatialgrid.z, np.pi)  # Extend the z grid to pi
        self.z_full_fine = np.linspace(self.z_full[0], self.z_full[-1], self.nz_extended)
        tortri = TorusCutTriangulation(common, z_res)
        self.triangulation = tortri.calc_triangulation()

    def generate_slice_attime(self, time):
        """ Specialised version of the slice version, interpolates onto the finer z grid"""
        slice3d_t = np.pad(super().generate_slice_attime(time), (0, 1), mode="wrap")
        pos_y = self._get_nearest_y_pos()
        posyrange = np.arange(0, self.cm.pnt.nky0*2 + 1)
        interp_data = np.empty((self.cm.pnt.nx0, self.nz_extended))

        for ix in range(len(self.cm.spatialgrid.x_a)):
            data_interp = interp.RectBivariateSpline(posyrange, self.z_full, slice3d_t[ix, ...])
            interp_data[ix, :] = data_interp(pos_y[ix, :], self.z_full_fine, grid=False)
        return interp_data

    def _get_nearest_y_pos(self):
        """ Calculate which y position we need for the toroidal slice """
        z = self.z_full_fine
        qprof = self.cm.pnt.q0*(
                1 + self.cm.pnt.shat/self.cm.pnt.x0*(self.cm.spatialgrid.x_a - self.cm.pnt.x0))
        # Get nearest y for given torangle \phi, y=C_y*(q(x)z - \phi)
        y_near = self.geom.Cy*(np.outer(qprof, z)*self.cm.pnt.sign_Ip_CW - self.torangle)
        ny = 2*self.cm.pnt.nky0
        dy = self.cm.pnt.rhostar*self.cm.pnt.ly/ny
        pos_y_near = (y_near/dy + ny/2)%ny
        return pos_y_near


class TorusCutTriangulation(object):
    """ Set up the triangulation of a torus cut based on GENE geometry outpout"""

    def __init__(self, common, z_res):
        """
        :param common: CommonData object of the run
        :param z_res: The number of points that z should be interpolated to
        """
        self.cm = common
        self.geom = self.geom = pydiag.utils.geom.Geometry(common)
        self.z_full = np.append(self.cm.spatialgrid.z, np.pi)  # Extend the z grid to pi
        self.z_full_fine = np.linspace(self.z_full[0], self.z_full[-1], z_res)
        self.triangulation = None

    def calc_triangulation(self, recalculate=False):
        """ Calculate the triangulation based on local/global

        :param recalculate: Force recalculation of the triangulation, otherwise a previous result is
        used
        """
        if not self.cm.y_local:
            raise RuntimeError("y-global is not supported in torus cut visualisation")
        if (not self.triangulation) or recalculate:
            if self.cm.x_local:
                return self._calc_triangulation_local()
            else:
                return self._calc_triangulation_xglobal()
        else:
            return self.triangulation

    def _calc_triangulation_xglobal(self):
        raise NotImplementedError("Torus viz for x-global is still work-in-progress")
        return None

    def _calc_triangulation_local(self):
        """ Calculate a triangulation based on R(x,z), Z(x,z) recovered from flux-tube geometry"""
        # These are the values at x=x0 extended so that z covers [-pi, pi]
        Rz = np.pad(self.geom.R, (0, 1), mode="wrap")
        Zz = np.pad(self.geom.Z, (0, 1), mode="wrap")
        # Interpolate for the x direction
        dxdR = np.pad(self.geom.dxdR/self.geom.gxx, (0, 1), mode="wrap")
        dxdZ = np.pad(self.geom.dxdZ/self.geom.gxx, (0, 1), mode="wrap")
        x_phys = (self.cm.spatialgrid.x_a - self.cm.pnt.x0)*self.cm.pnt.Lref
        R_pos = interp.RectBivariateSpline(x_phys, self.z_full, Rz + np.outer(x_phys, dxdR))(x_phys,
                                                                                             self.z_full_fine)
        Z_pos = interp.RectBivariateSpline(x_phys, self.z_full, Zz + np.outer(x_phys, dxdZ))(x_phys,
                                                                                             self.z_full_fine)
        # Calculate a triangulation based on R(x,z) and Z(x,z)
        RZtri = mtri.Triangulation(R_pos.flatten(), Z_pos.flatten())
        # Mask central hole
        centretri = mtri.Triangulation(R_pos[0, :], Z_pos[0, :])
        centermask = np.full(RZtri.triangles.shape[0], False)
        for ct in centretri.triangles:
            start = np.where(ct.mean() == RZtri.triangles.mean(axis=1))
            centermask[start] = True
        RZtri.set_mask(centermask)
        self.triangulation = RZtri
        return self.triangulation


class PlotTorusCut(pydiag.diagplots.baseplot.Plotting):
    """ Class handling the plotting of the triangulated torus cuts"""

    def __init__(self, toruscutlist):
        super().__init__()
        self.toruscutlist = toruscutlist
        self.anims = []  # The animation objects need to be stored for the animations to run
        self.figlist = []
        self.ax = None

    def createfigs(self, animate=False, save=False):
        for torus in self.toruscutlist:
            if not animate:
                torus.generate_timeaverage()
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_aspect('equal', 'box')
                ax.tricontourf(torus.triangulation, torus.timeaverage.flatten(), 50,
                               cmap=self.cmap_unidirect)
                ax.set_xlabel(r"$R[m]$")
                ax.set_ylabel(r"$Z[m]$")
                fig.tight_layout()
            else:
                torus.generate_timeseries()
                maxval = np.max(torus.dataarray)
                minval = np.min(torus.dataarray)
                self.figlist.append(plt.figure())
                self.ax = self.figlist[-1].add_subplot(111)
                self.ax.set_aspect('equal', 'box')
                self.ax.tricontourf(torus.triangulation, torus.dataarray[0].flatten(), 50,
                                    cmap=self.cmap_unidirect, vmax=maxval, vmin=minval, zorder=-20)

                def on_click(event):
                    if event.dblclick:
                        for ani in self.anims:
                            if ani.running:
                                ani.event_source.stop()
                            else:
                                ani.event_source.start()
                            ani.running ^= True

                def animate(t):
                    self.ax.clear()
                    self.ax.set_xlabel(r"$R[m]$")
                    self.ax.set_ylabel(r"$Z[m]$")
                    self.ax.set_rasterization_zorder(z=-10)
                    self.ax.set_title(r"t={:.3f}".format(torus.timearray[t]))
                    cma = self.ax.tricontourf(torus.triangulation, torus.dataarray[t].flatten(), 50,
                                              cmap=self.cmap_unidirect, vmax=maxval, vmin=minval,
                                              zorder=-20)
                    return cma

                self.anims.append(animation.FuncAnimation(self.figlist[-1], animate,
                                                          range(len(torus.timearray)), blit=False))
                for anim in self.anims:
                    anim.running = True
                for fig in self.figlist:
                    fig.canvas.mpl_connect('button_press_event', on_click)
                if save:
                    # Keep the plotting from freezing by doing the saving in a subprocess
                    def save_anim():
                        print("Saving animation in the background, this can take a long while")
                        self.anims[-1].save(
                                "animatedtorus_{}_{}{}.mp4".format(torus.quantity, torus.species,
                                                                   torus.cm.fileextension), dpi=300)
                        print("Finished saving animation")

                    p = multiprocessing.Process(target=save_anim)
                    p.daemon = True
                    p.start()
