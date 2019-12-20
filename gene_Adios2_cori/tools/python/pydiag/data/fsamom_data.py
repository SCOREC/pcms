"""fsamom_data.py

Contains the infrastructure to handle fsamom_spec_fext files

"""
import numpy as np
from bisect import bisect_left
from pydiag.data.base_file import TimeSeries
from pydiag.utils.ParIO import Parameters


class FSAmomFile(TimeSeries):
    """Class to handle a specific fsamom file

    :param filename: Filename, typically fsamom[n,p,e]_[species]_fext
    :param common: CommonData object of the run
    """

    def __init__(self, filename, common):
        super().__init__(filename, common)
        self.fsacolumns = 13  # Number of columns in fsamom file
        self.dataarray = []
        self.momtitle = []
        self.generate_timeseries()

    def generate_timeseries(self):
        """ Read flux surface averaged moments """
        blocks = []
        self.timearray = []
        print('Reading {}\n'.format(self.objectname))
        try:
            with open(self.objectname) as fsafile:
                next(fsafile)  # skip header
                for line in fsafile:
                    if not line or line.startswith('\n'):
                        continue
                    if line.startswith("#   x/a"):
                        continue
                    if line.startswith('#'):
                        self.timearray.append(float(line.split()[1]))
                        blocks.append([])
                    else:
                        blocks[-1].append(line)
        except IOError:
            print("IOError: probably fsa file does not exist: {}".format(self.objectname))
            raise
        self.timearray = np.array(self.timearray)
        self.check_times()
        pos = self.calc_positions()
        self._addfsadata(blocks, pos)
        self.timearray = self.timearray[pos]

    def _addfsadata(self, blocks, pos):
        for entry in pos:
            self.dataarray.append(np.empty((self.cm.pnt.nx0, self.fsacolumns)))
            for i in range(0, self.cm.pnt.nx0):
                line = blocks[entry][i]
                for var in range(self.fsacolumns):
                    try:
                        self.dataarray[-1][i, var] = float(line.split()[var])
                    except ValueError:
                        self.dataarray[-1][i, var] = 0.0


class FSAmomData(object):
    """Class to handle the flux surface averaged moment diag of a specific run

    Version that handles the new version of the diagnostic with three output files
    """

    def __init__(self, common):
        self.isdatapresent = False
        self.cm = common
        self.fsaspec = []

    def getfsadata(self):
        """Read the data from the fsa_moments diagnostic"""
        try:
            if self.cm.pnt.istep_fsa_moments == 0:
                print("No fsamom diagnostic included in GENE run")
                return
        except AttributeError:
            print("No fsamom diagnostic included in GENE run")
            return
        for n in range(0, self.cm.pnt.n_spec):
            filename = 'fsamomn_{}{}'.format(self.cm.specnames[n], self.cm.fileextension)
            self.fsaspec.append(FSAmomFile(filename, self.cm))
            self.fsaspec[-1].momtitle = r"Density"
            filename = 'fsamomp_{}{}'.format(self.cm.specnames[n], self.cm.fileextension)
            self.fsaspec.append(FSAmomFile(filename, self.cm))
            self.fsaspec[-1].momtitle = r"$v_\parallel$"
            filename = 'fsamome_{}{}'.format(self.cm.specnames[n], self.cm.fileextension)
            self.fsaspec.append(FSAmomFile(filename, self.cm))
            self.fsaspec[-1].momtitle = r"$v^2$"
        self.isdatapresent = True


class FSAmomData_legacy(object):
    """Class to handle the flux surface averaged moment diag of a specific run

    Version that handles the old version of the diagnostic with one output files
    """

    def __init__(self, common):
        self.isdatapresent = False
        self.cm = common
        self.fsaspec = []

    def getfsadata(self):
        """Read the data from the fsa_moments diagnostic"""
        try:
            if self.cm.pnt.istep_fsa_moments == 0:
                print("No fsamom diagnostic included in GENE run")
                return
        except AttributeError:
            print("No fsamom diagnostic included in GENE run")
            return
        for n in range(0, self.cm.pnt.n_spec):
            filename = 'fsamom_{}{}'.format(self.cm.specnames[n], self.cm.fileextension)
            self.fsaspec.append(FSAmomFile(filename, self.cm))
            self.fsaspec[-1].momtitle = r"Density, $v_\parallel$ or $v^2$"
        self.isdatapresent = True
