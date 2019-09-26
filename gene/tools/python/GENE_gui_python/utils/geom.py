"""Module containing the treatment of GENE's geometry output"""
import numpy as np
import h5py


class Geometry:
    """ Class to handle geometry input from GENE runs

    """

    def __init__(self, run):
        self.geomtype = run.pars['magn_geometry']
        if run.is_h5:
            self.getgeom_h5(run)
        elif run.is_adios:
            raise NotImplementedError("ADIOS not supported")
        else:
            self.getgeom_std(run)

    def getgeom_std(self, run):
        self.cm = run
        geom = self.getgeom(run)
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
        self.dxdR = geom[13]
        self.dxdZ = geom[14]

        if not run.x_local:
            self.Cy = geom[15, 0]
            self.Cxy = geom[16, 0]
            self.q = geom[17, 0]
            self.dpdx_pm_arr = geom[18, 0]
            self.jaco3d = np.broadcast_to(self.jacobian[:, np.newaxis, :],
                                          (run.pnt.nz0, run.pnt.nky0, run.pnt.nx0))
        elif not run.y_local:
            self.jaco3d = np.broadcast_to(self.jacobian[:, :, np.newaxis],
                                          (run.pnt.nz0, run.pnt.nky0, run.pnt.nx0))

    @staticmethod
    def untangle_1d(arrlist, start, nxorny):
        """ Process a geometry file section that contains a 1d (x or y) field

        :param arrlist: The entire file as a list of lines
        :param start: The starting index of the quantity in arrlist
        :param nxorny: nx0 or nky0 depending on the type of simulation
        :returns: A 1d numpy array
        """
        arr_1d = np.zeros(nxorny)
        i = 0
        for line in arrlist[start:-1]:
            for j in range(len(line.split())):
                arr_1d[i + j] = line.split()[j]
            i += len(line.split())
            if i >= nxorny:
                break
        return arr_1d

    @staticmethod
    def untangle_2d(arrlist, start, nxorny, nz):
        """ Process a geometry file section that contains a 2d (x or y, z) field

        :param arrlist: The entire file as a list of lines
        :param start: The starting index of the quantity in arrlist
        :param nxorny: nx0 or nky0 depending on the type of simulation
        :param nz: nz0, number of z grid points
        :returns: A 2d numpy array
        """
        arr_2d = np.zeros((nz, nxorny))
        ik = 0
        for line in arrlist[start:-1]:
            arr = line.split()[:]
            lb = ik%nxorny
            ub = min(lb + 16, nxorny)
            k = int(ik/nxorny)
            dim = ub - lb
            arr_2d[k, lb:ub] = arr[0:dim]
            if dim < 16:
                arr_2d[k + 1, 0:16 - dim] = arr[dim:]
            ik += 16
            if ik >= nz*nxorny:
                break
        return arr_2d

    def getgeom(self, run):
        """ Returns the geometry from a non-hdf5 file """
        local = run.x_local and run.y_local
        with open(FileName(run.folder, self.geomtype.strip("'"), run.fileextension),
                  "r") as geomfile:
            if local:
                geom = np.empty((16, run.pnt.nz0), dtype=np.float64)
                k = 0
                for line in geomfile:
                    if len(line.split()) == 16:
                        geom[:, k] = line.split()[:]
                        k += 1
                    elif line.startswith('Cy'):
                        self.Cy = float(line.split()[-1])
                    elif line.startswith('Cxy'):
                        self.Cxy = float(line.split()[-1])
            elif run.y_local:  # x-global
                geom = self.getgeom_glob(geomfile, run.pnt.nz0, run.pnt.nx0)
            elif run.x_local:  # y-global
                geom = self.getgeom_glob(geomfile, run.pnt.nz0, run.pnt.nky0)
            else:
                raise NotImplementedError("xy global not supported")
        return geom

    def getgeom_glob(self, geomfile, nz0, nxorny):
        """ Subroutine for geometry files from global runs """
        NUMFIELDS = 19
        geom = np.zeros((NUMFIELDS, nz0, nxorny), dtype=np.float64)

        geomdict = {'gxx': 0, 'gxy': 1, 'gxz': 2, 'gyy': 3, 'gyz': 4, 'gzz': 5, 'Bfield': 6,
                    'dBdx': 7, 'dBdy': 8, 'dBdz': 9, 'jacobian': 10, 'geo_R': 11, 'geo_Z': 12,
                    'geo_c1': 13, 'geo_c2': 14, 'C_y': 15, 'C_xy': 16, 'q': 17, 'dpdx_pm_arr': 18}
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
                    pos_start[num] = linnum + 1
                except KeyError:
                    try:  # Test if it is a single number (can occur for 1d arrays)
                        float(line.strip())
                    except ValueError:
                        raise RuntimeError("Unknown entry name in geometry file")
        for num in range(NUMFIELDS):
            if pos_start[num] == 0:  # This should only occur if a field does not exist
                continue
            if num in [15, 16, 17, 18]:  # We have a 1d field following
                geom[num] = self.untangle_1d(geomlist, pos_start[num], nxorny)
            else:
                geom[num] = self.untangle_2d(geomlist, pos_start[num], nxorny, nz0)

        return geom

    def getgeom_h5(self, run):
        """ Returns the geometry from a non-hdf5 file """

        geomfile = FileName(run.folder, self.geomtype.strip("'"), run.fileextension)
        self.cm = run
        geom = h5py.File(geomfile, 'r')

        self.gxx = geom.get('/metric/g^xx').value
        self.gxy = geom.get('/metric/g^xy').value
        self.gxz = geom.get('/metric/g^xz').value
        self.gyy = geom.get('/metric/g^yy').value
        self.gyz = geom.get('/metric/g^yz').value
        self.gzz = geom.get('/metric/g^zz').value
        self.Cy = geom.get('/metric/C_y').value
        self.Cxy = geom.get('/metric/C_xy').value

        self.Bfield = geom.get('/Bfield_terms/Bfield').value
        self.dBdx = geom.get('/Bfield_terms/dBdx').value
        self.dBdy = geom.get('/Bfield_terms/dBdy').value
        self.dBdz = geom.get('/Bfield_terms/dBdz').value
        self.jacobian = geom.get('/Bfield_terms/Jacobian').value

        if not (not run.x_local and not run.y_local):
            self.R = geom.get('/shape/R').value
            self.Z = geom.get('/shape/Z').value
            self.dxdR = geom.get('/shape/dxdR').value
            self.dxdZ = geom.get('/shape/dxdZ').value

        try:
            self.q = geom.get('/profile/q_prof').value
        except:
            self.q = run.pnt.q0

        try:
            self.dpdx_arr = geom.get('/profile/dpdx_pm_arr').value
        except:
            self.dpdx_arr = run.pnt.dpdx_pm

        geom.close()


def FileName(folder, prefix, ext):
    return check_str(folder) + prefix + ext


def check_str(folder):
    if folder[-1:] != '/':
        folder = folder + '/'
    return folder
