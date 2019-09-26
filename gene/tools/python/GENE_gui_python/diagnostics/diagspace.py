#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

    def __init__(self, rungrid, x_fourier, y_fourier, z_fourier=False, xrange=(None, None),
                 yrange=(None, None), zrange=(None, None), xavg=False, yavg=False, zavg=False):
        """ Standard constructor for DiagSpace

        :param rungrid: SpatialGrid object, supplies the numerical grid of the GENE run
        :param x_fourier: If true, operate in kx space, Fourier transform data if necessary
        :param y_fourier: If true, operate in ky space, Fourier transform data if necessary
        :param z_fourier: If true, operate in kz space, Fourier transform data if necessary
        :param xrange: x or kx index range as a tuple (i_min, i_max). Note that kx refers to the
        positive kx modes here.
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
        self.z_fourier = z_fourier
        if x_fourier:
            self.xslice = self._create_slice(xrange, rungrid.kx_pos)
        else:
            self.xslice = self._create_slice(xrange, rungrid.x)
        if y_fourier:
            self.yslice = self._create_slice(yrange, rungrid.ky)
        else:
            self.yslice = self._create_slice(yrange, rungrid.y)
        if z_fourier:
            self.zslice = self._create_slice(zrange, rungrid.kz)
        else:
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
        self.rungrid = rungrid

    @staticmethod
    def _create_slice(coordrange, grid_direction):
        if grid_direction[coordrange[0]:coordrange[1]].size == 0:
            raise IndexError(
                    "The given range {} is empty or outside of the grid ({}:{}) of the run".format(
                            coordrange, 0, grid_direction.size - 1))
        return slice(coordrange[0], coordrange[1])
