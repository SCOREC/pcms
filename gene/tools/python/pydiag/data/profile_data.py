""" profile_data.py

Contains the infrastructure to handle profile_spec_fileextension files
Also includes nustar and the Chang-Hinton neoclassical heat flux prediction
"""

import math
import numpy as np

from pydiag.data.base_file import TimeSeries
from pydiag.utils.geom import Geometry
from pydiag.utils.nc_predictions import chang_hinton


class ProfileFile(TimeSeries):
    """Class to handle a profile diagnostic file

    :param filename: Typically profile_{spec}{fileextension}
    :param common: CommonData object of the run
    """

    def __init__(self, filename, common):
        super().__init__(filename, common)
        self.dataarray = {"Ts": [], "ns": [], "omts": [], "omns": [], "Gammanc": [], "Gammaturb": [],
                         "Qturb": [], "Qnc": [], "Pinc": [], "Piturb": [], "jbs": []}
        # Collisionality, i.e. ratio of \nu and bounce frequency
        self.nustar = []

    def generate_timeseries(self):
        """ Read density, temperature, gradients and fluxes """
        self.timearray, blocks = self._parse_profile_file()
        self.timearray = np.array(self.timearray)
        self.check_times()
        #  For a single time, find element closest to given input time
        pos = self.calc_positions()
        self._addprofdata(blocks, pos)
        self.timearray = self.timearray[pos]

    def _parse_profile_file(self):
        blocks = []
        timearray = []
        print('Reading  {}\n'.format(self.objectname))
        try:
            with open(self.objectname) as prfile:
                next(prfile)  # skip header
                for line in prfile:
                    if not line or line.startswith('\n'):
                        continue
                    if line.startswith('#'):
                        timearray.append(float(line.split()[1]))
                        blocks.append([])
                    else:
                        blocks[-1].append(line)
        except IOError:
            print("Probably profile file does not exist: {}".format(self.objectname))
            raise
        return timearray, blocks

    def _addprofdata(self, blocks, pos):
        """ Fill vectors with profile data """
        for entry in pos:
            self.dataarray["Ts"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["ns"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["omts"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["omns"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Gammaturb"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Qturb"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Piturb"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Gammanc"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Qnc"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["Pinc"].append(np.zeros(self.cm.pnt.nx0))
            self.dataarray["jbs"].append(np.zeros(self.cm.pnt.nx0))
            for i in range(0, self.cm.pnt.nx0):
                linearr = (blocks[entry][i]).split()
                self.dataarray["Ts"][-1][i] = linearr[2]
                self.dataarray["ns"][-1][i] = linearr[3]
                self.dataarray["omts"][-1][i] = linearr[4]
                self.dataarray["omns"][-1][i] = linearr[5]
                self.dataarray["Gammaturb"][-1][i] = linearr[6]
                self.dataarray["Qturb"][-1][i] = linearr[7]
                self.dataarray["Piturb"][-1][i] = linearr[8]
                self.dataarray["Gammanc"][-1][i] = linearr[9]
                self.dataarray["Qnc"][-1][i] = linearr[10]
                self.dataarray["Pinc"][-1][i] = linearr[11]
                self.dataarray["jbs"][-1][i] = linearr[12]
            # Correct for convective energy flux
            self.dataarray["Qnc"][-1] -= 2.5*self.dataarray["Gammanc"][-1]*self.dataarray["Ts"][-1]
            self.dataarray["Qturb"][-1] -= 2.5*self.dataarray["Gammaturb"][-1]*self.dataarray["Ts"][-1]


class ProfileData(object):
    """Class handling the profile diagnostic files of a run (i.e. one per species)

    :param common: CommonData object of the run
    """

    def __init__(self, common, rundatafiles):
        """ Read parameter file and create empty arrays for profile data """
        self.isdatapresent = False
        self.cm = common
        self.rundatafiles = rundatafiles
        self.geom = Geometry(common)
        self.q = np.array(self.geom.q)
        self.prspec = []
        self.xs = np.empty(self.cm.pnt.nx0)
        self.T0s = np.empty((self.cm.pnt.nx0, self.cm.pnt.n_spec))
        self.n0s = np.empty((self.cm.pnt.nx0, self.cm.pnt.n_spec))
        self.omt0s = np.empty((self.cm.pnt.nx0, self.cm.pnt.n_spec))
        self.omn0s = np.empty((self.cm.pnt.nx0, self.cm.pnt.n_spec))

        self.Afs = []  # Flux surface area
        # Chang-Hinton prediction for main ion species
        self.Qcharr = []

    def calc_afs(self):
        """ Calculate area of flux surfaces"""
        c = self.cm.pnt.ly*self.cm.pnt.rhostar*self.cm.pnt.minor_r*2.0*np.pi*self.cm.pnt\
            .n0_global/self.cm.pnt.nz0
        self.Afs = c*np.sum(self.geom.jacobian*np.sqrt(self.geom.gxx), axis=0)/self.cm.pnt.major_R

    def _read_zero_prof(self):
        """Get the initial profiles from profiles_spec_fileextension"""
        for n in range(self.cm.pnt.n_spec):
            prof0file = open('profiles_{}{}'.format(self.cm.specnames[n], self.cm.fileextension))
            lines = prof0file.readlines()
            for i in range(0, self.cm.pnt.nx0):
                l = 2 + i
                self.xs[i] = lines[l].split()[0]
                self.T0s[i, n] = float(lines[l].split()[2]) / self.cm.pnt.Tref
                self.n0s[i, n] = float(lines[l].split()[3]) / self.cm.pnt.nref
                self.omt0s[i, n] = lines[l].split()[4]
                self.omn0s[i, n] = lines[l].split()[5]
            prof0file.close()

    def get_profiles(self):
        """ Organize the profile data, calls most other functions"""
        self._read_zero_prof()
        for n in range(0, self.cm.pnt.n_spec):
            prs = self.rundatafiles.get_fileobject('profile_{}'.format(self.cm.specnames[n]))
            prs.generate_timeseries()
            self.prspec.append(prs)
        self.calc_nustar()
        self.Qcharr = chang_hinton(self, 0)  # Only main ion species
        self.calc_afs()
        self.isdatapresent = True

    def calc_nustar(self):
        """ Calculate the collisionality nustar for each species. """
        try:
            coll = self.cm.pnt.coll
        except AttributeError:
            coll = 0
        epsilon = np.array(self.xs * self.cm.pnt.minor_r / self.cm.pnt.major_R)
        for prs in self.prspec:
            prs.nustar = []
            for n, T in zip(prs.dataarray["ns"], prs.dataarray["Ts"]):
                prs.nustar.append(np.array(8.0 / 3.0 / math.sqrt(math.pi) *
                                           self.cm.pnt.major_R * coll / epsilon ** (3.0 / 2) *
                                           self.q * n / T ** 2))

    def generate_timetrace(self, xpos, quantity="Qturb"):
        """ Generate a 1d-timetrace of a profile quantity.

        :param xpos: Radial position in x/a, the nearest x gridpoint will be used.
        :param quantity: Determines the quantity which will be used, default turbulent heat flux
        """
        xind = (np.abs(self.xs - xpos)).argmin()  # find index of position closest to xpos
        timetr = []
        for n, prs in enumerate(self.prspec):
            for step in prs.dataarray[quantity]:
                timetr.append(step[xind])
            timetrfile = open("ttrace_{}{}_{}{}".format(xpos, quantity, self.cm.specnames[n],
                                                        self.cm.fileextension), mode='w')
            timetrfile.write("# time  " + quantity + "\n")
            for t, q in zip(prs.timearray, timetr):
                timetrfile.write(str(t) + "  " + str(q) + "\n")
            timetrfile.close()


def glueprofiledata(prdlist):
    """ Take profdata objects with matching time intervals and combine them into a larger one
    :param prdlist: List of ProfileData objects to combine
    :returns: Combined Profdata object

    """
    # TODO: sanity checks for the parameters of prdlist entries
    # Nothing to combine, nothing to do
    if len(prdlist) == 1:
        return prdlist[0]
    # Use first ProfileData as basis for the construction of the glued object
    result = prdlist.pop(0)
    for prd in prdlist:
        for n in range(result.cm.pnt.n_spec):
            # When times overlap, give following ProfileData preference
            # If the timefield has been converted to a np array before we need to undo it
            result.prspec[n].timearray = result.prspec[n].timearray.tolist()
            prd.prspec[n].timearray = prd.prspec[n].timearray.tolist()

            while result.prspec[n].timearray[-1] >= prd.prspec[n].timearray[0]:
                del result.prspec[n].timearray[-1]
                for key in result.prspec[n].dataarray:
                    del result.prspec[n].dataarray[key][-1]
            result.prspec[n].timearray = np.concatenate((result.prspec[n].timearray,
                                                         prd.prspec[n].timearray))
            for key in result.prspec[n].dataarray:
                result.prspec[n].dataarray[key] += prd.prspec[n].dataarray[key]
    result.cm.endtime = prdlist[-1].cm.endtime
    result.calc_nustar()
    result.Qcharr = chang_hinton(result, 0)
    del prdlist
    return result
