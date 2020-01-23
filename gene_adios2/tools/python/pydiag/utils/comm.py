import numpy as np

from pydiag.utils.ParIO import Parameters


class CommonData(object):
    """ Class containing basic information shared between diagnostics

    :param fileextension: File extension of the run to analyse
    :param starttime: Start of the analysis time window
    :param endtime: End of the analysis time window

    Special choices for starttime and endtime:
    -1: use first time in profile file
    -2: use last time in profile file
    """

    def __init__(self, fileextension, starttime, endtime):
        self.starttime = starttime
        self.endtime = endtime
        self.fileextension = fileextension
        params = Parameters()
        params.Read_Pars('parameters{}'.format(fileextension))
        self.pars = params.pardict
        self.pnt = params.asnamedtuple()
        # Determine if we are in one of the global versions
        self.x_local = self.pars.get('x_local', True)  # By default, x_local = T is not written
        self.y_local = self.pars.get('y_local', True)  # By default, y_local = T is not written
        self.nonlinear = self.pars.get('nonlinear')
        self.specnames = [str(self.pars['name{:1d}'.format(nsp + 1)])[1:-1] for nsp in
                          range(self.pnt.n_spec)]
        self.spatialgrid = SpatialGrid(self.pnt, self.x_local, self.y_local)


class SpatialGrid(object):
    """ The numerical grid based on the parameters of a GENE run

    :param pnt: NamedTuple of the run parameters
    :param x_local: Is the run local in radial direction?
    :param y_local: Is the run local in "binormal" direction?
    """

    def __init__(self, pnt, x_local, y_local):
        self._calc_kxgrid(x_local, y_local, pnt.lx, pnt.nx0)
        self._calc_xgrid(x_local, pnt.lx, pnt.nx0, pnt.rhostar, pnt.x0)
        self._calc_kygrid(y_local, pnt.kymin, pnt.nky0)
        self._calc_ygrid(y_local, pnt.kymin, pnt.nky0)
        self.z = np.array([-np.pi + 2*np.pi*k/pnt.nz0 for k in range(pnt.nz0)], dtype=np.float64)

    def _calc_xgrid(self, x_local, lx, nx0, rhostar, x0):
        if x_local:
            xgrid = np.fft.fftfreq(nx0, d=1.0/lx)
            self.x = np.fft.fftshift(xgrid)
        else:
            self.x = np.array([-lx/2. + lx/(nx0 - 1)*i for i in range(0, nx0)])
            self.x_a = self.x*rhostar + x0

    def _calc_kxgrid(self, x_local, y_local, lx, nx0):
        self.kxmin = 2*np.pi/lx
        if x_local:
            kxgridp = [self.kxmin*i for i in range(int(nx0/2 + 1))]
            kxgridn = [self.kxmin*i for i in range(int(-nx0/2 + 1), 0)]
            if y_local:
                self.kx_fftorder = np.array(kxgridp + kxgridn, dtype=np.float64)
                # ordered kx array
                self.kx = np.array(kxgridn + kxgridp[0:-1], dtype=np.float64)
            else:
                self.kx = np.array([self.kxmin*i for i in range(nx0)])
        else:
            if y_local:
                self.kx_fftorder = 2*np.pi*np.fft.fftfreq(nx0, d=lx/nx0)
                self.kx = np.fft.ifftshift(self.kx_fftorder)
            else:
                self.kx_fftorder = 2*np.pi*np.fft.fftfreq(nx0, d=lx/2/nx0)
                self.kx = np.fft.ifftshift(self.kx_fftorder)

    def _calc_ygrid(self, y_local, kymin, nky0):
        ly = 2*np.pi/kymin
        if y_local:
            self.y = np.fft.fftshift(np.fft.fftfreq(2*nky0, d=1/ly))
        else:
            self.y = np.array([-ly/2. + ly/nky0*i for i in range(0, nky0)])

    def _calc_kygrid(self, y_local, kymin, nky0):
        if y_local:
            self.ky = np.array([kymin*j for j in range(nky0)], dtype=np.float64)
        else:
            ly = 2*np.pi/kymin
            self.ky = 2*np.pi*np.fft.ifftshift(np.fft.fftfreq(nky0, d=ly/nky0))


class DiagSpace(object):
    """ Class defining the space that should be used in a diagnostic

    i.e. which intervals in x, y, z, Fourier or configuration space in x,y
    Velocity space is not handled by this class
    The ranges are given as tuples of the starting and stopping index (which is not included)
    Negative indices can be used in the usual Python way (-1 indexing the last element).
    Be aware that while [0:] includes the last element [0:-1] does not!
    Refer to https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html for the details
    how indexing numpy arrays works.
    A binary flag indicates if an average is desired in that direction.
    """

    def __init__(self, rungrid, x_fourier, y_fourier, xrange=(None, None), yrange=(None, None),
                 zrange=(None, None), xavg=False, yavg=False, zavg=False):
        """ Standard constructor for DiagSpace

        :param rungrid: SpatialGrid object, supplies the numerical grid of the GENE run
        :param x_fourier: If true, operate in kx space, Fourier transform data if necessary
        :param y_fourier: If true, operate in ky space, Fourier transform data if necessary
        :param xrange: x or kx index range as a tuple (i_min, i_max)
        :param yrange: y or ky index range as a tuple (j_min, j_max)
        :param zrange: z index range as a tuple (k_min, k_max)

        Note that None is also a possible choice for the min and max values due to [:]
        giving the entire array.
        """
        self.xavg = xavg
        self.yavg = yavg
        self.zavg = zavg
        self.x_fourier = x_fourier
        self.y_fourier = y_fourier
        if x_fourier:
            self.xslice = self._create_slice(xrange, rungrid.kx)
        else:
            self.xslice = self._create_slice(xrange, rungrid.x)
        if y_fourier:
            self.yslice = self._create_slice(yrange, rungrid.ky)
        else:
            self.yslice = self._create_slice(yrange, rungrid.y)
        self.zslice = self._create_slice(zrange, rungrid.z)
        self.xyslice = (self.xslice, self.yslice)
        self.xzslice = (self.xslice, self.zslice)
        self.yzslice = (self.yslice, self.zslice)
        if not (xavg or yavg or zavg):
            self.diagslice = (self.xslice, self.yslice, self.zslice)
        elif xavg and not yavg and not zavg:
            self.diagslice = self.yzslice
        elif xavg and yavg and not zavg:
            self.diagslice = self.zslice
        elif xavg and not yavg and zavg:
            self.diagslice = self.yslice
        elif not xavg and not yavg and zavg:
            self.diagslice = self.xyslice
        elif not xavg and yavg and not zavg:
            self.diagslice = self.xzslice
        elif not xavg and yavg and zavg:
            self.diagslice = self.xslice
        else:
            self.diagslice = None

    @staticmethod
    def _create_slice(coordrange, grid_direction):
        if grid_direction[coordrange[0]:coordrange[1]].size == 0:
            raise IndexError(
                    "The given range {} is empty or outside of the grid ({}:{}) of the run".format(
                        coordrange, 0, grid_direction.size-1))
        return slice(coordrange[0], coordrange[1])
