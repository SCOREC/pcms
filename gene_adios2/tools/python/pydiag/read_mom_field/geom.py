"""Module containing the treatment of GENE's geometry output"""
# pylint: disable=E1101
import numpy as np
from os.path import join


class geometry:
    """ Class to handle geometry input from GENE runs """
    def __init__(self, cm, p):
        global common, pars
        common = cm
        pars = p

        self.geomtype = pars.pardict['magn_geometry']
        geom = self.getgeom()
        self.gxx = geom[0]
        self.gxy = geom[1]
        self.gxz = geom[2]
        self.gyy = geom[3]
        self.gyz = geom[4]
        self.gzz = geom[5]
        self.Bfield = geom[6]
        self.dBdx = geom[7]
        self.dBdy = geom[8]
        self.dBdz = geom[9]
        self.jacobian = geom[10]
        self.R = geom[11]
        self.Z = geom[12]
        self.drhodr = geom[13]
        self.drhodz = geom[14]
        self.Cy = geom[15, 0]
        if not common.x_local:
            self.Cxy = geom[16, 0]
            self.q = geom[17, 0]
            self.dpdx_pm_arr = geom[18, 0]

    @staticmethod
    def untangle_1d(arrlist, start, nxorny):
        """ Process a geometry file section that contains a 1d (x or y) field """
        arr_1d = np.zeros(nxorny)
        i = 0
        for line in arrlist[start:-1]:
            for j in range(len(line.split())):
                arr_1d[i+j] = line.split()[j]
            i += len(line.split())
            if i >= nxorny:
                break
        return arr_1d

    @staticmethod
    def untangle_2d(arrlist, start, nxorny, nz):
        """ Process a geometry file section that contains a 2d (x or y, z) field """
        arr_2d = np.zeros((nz, nxorny))
        ik = 0
        for line in arrlist[start:-1]:
            arr = line.split()[:]
            lb = ik%nxorny
            ub = min(lb+16, nxorny)
            k = int(ik/nxorny)
            dim = ub-lb
            arr_2d[k, lb:ub] = arr[0:dim]
            if dim < 16:
                arr_2d[k+1, 0:16-dim] = arr[dim:]
            ik += 16
            if ik >= nz*nxorny:
                break
        return arr_2d

    def getgeom(self):
        """ Returns the geometry from a non-hdf5 file """
        xlocal = False
        ylocal = False
        try:
            if pars.pardict['x_local'] in ['.t.', '.T.', 't']:
                xlocal = True
        except KeyError:
            xlocal = True
        try:
            if pars.pardict['y_local'] in ['.t.', '.T.', 't']:
                ylocal = True
        except KeyError:
            ylocal = True
        local = xlocal and ylocal
        with open(join(common.rundir, self.geomtype.strip("'"))+common.fileext, "r") as geomfile:
            if local:
                geom = np.empty((16, common.nz), dtype=np.float64)
                k = 0
                for line in geomfile:
                    if len(line.split()) == 16:
                        geom[:, k] = line.split()[:]
                        k += 1
                    elif line.startswith('Cxy'):
                        self.Cxy = float(line.split()[-1])
            elif ylocal:  # x-global
                geom = self.getgeom_glob(geomfile, common.nx)
            elif xlocal:  # y-global
                geom = self.getgeom_glob(geomfile, common.ny)
            else:
                raise NotImplementedError("xy global not supported")
        return geom

    def getgeom_glob(self, geomfile, nxorny):
        """ Subroutine for geometry files from global runs """
        NUMFIELDS = 19
        geom = np.zeros((NUMFIELDS, common.nz, nxorny), dtype=np.float64)

        geomdict = {'gxx': 0, 'gxy': 1, 'gxz': 2, 'gyy': 3, 'gyz': 4, 'gzz': 5,
                    'Bfield': 6, 'dBdx': 7, 'dBdy': 8, 'dBdz': 9, 'jacobian': 10,
                    'geo_R': 11, 'geo_Z': 12, 'geo_c1': 13, 'geo_c2': 14, 'C_y': 15,
                    'C_xy': 16, 'q': 17, 'dpdx_pm_arr': 18}
        pos_start = np.zeros(NUMFIELDS, dtype=int)  # where a field starts
        num = 0
        parheader = True
        geomlist = geomfile.readlines()
        for linnum, line in enumerate(geomlist):
            if parheader:  # ignore the top part of the file (copy of the geometry namelist)
                if line.startswith(r'/'):
                    parheader = False
                continue
            if len(line.split()) == 1:  # the variable names
                try:
                    num = geomdict[line.strip()]
                    pos_start[num] = linnum+1
                except KeyError:
                    try:  # Test if it is a single number (can occur for 1d arrays)
                        float(line.strip())
                    except ValueError:
                        raise RuntimeError("Unknown entry name in geometry file")
        for num in range(NUMFIELDS):
            if pos_start[num] == 0:     # This should only occur if a field does not exist
                continue
            if num in [15, 16, 17, 18]:  # We have a 1d field following
                geom[num] = self.untangle_1d(geomlist, pos_start[num], nxorny)
            else:
                geom[num] = self.untangle_2d(geomlist, pos_start[num], nxorny, common.nz)
        return geom
