"""
Module handling nrg files

"""
import csv
import numpy as np

from pydiag.data.base_file import TimeSeries


class NrgFile(TimeSeries):
    """NrgFile: Class to read a nrg file from GENE

    :param filename: Typically nrg.dat
    :param common: CommonData object of the run
    """
    def __init__(self, filename, common):
        super().__init__(filename, common)
        self.dataarray = []
        self.n_col = common.pnt.nrgcols
        self.n_spec = common.pnt.n_spec

    def generate_timeseries(self):
        """ Fill the Nrgdata object with data """
        n_col = self.n_col
        self.timearray = []
        self.dataarray = []
        try:
            with open(self.objectname) as nrgfile:
                csvnrg = csv.reader(nrgfile, delimiter=' ', skipinitialspace=True)
                for line in csvnrg:
                    if len(line) == 0:
                        continue
                    if len(line) == 1:
                        self.timearray.append(float(line[0]))
                        self.dataarray.append([[] for _ in range(self.n_spec)])
                        ispec = 0  # Reset index for the species at the current time step
                    elif len(line) == n_col:
                        self.dataarray[-1][ispec] = line
                        ispec += 1
                    else:
                        raise IOError("Incorrect number of columns")
        except IOError:
            raise IOError("nrg file does not exist or has"
                          " wrong number of columns: {}".format(self.objectname))
        self.check_times()
        pos = self.calc_positions()
        self.dataarray = np.array(self.dataarray).astype(float, copy=False)
        # Reduce the nrgcols to only the required time frame
        self.dataarray = self.dataarray[pos, ...]
        self.timearray = np.array(self.timearray)[pos]


def gluenrgdata(nrglist):
    """ Function to combine a list of nrg files from a continuation run.

    :param nrglist: the list of NrgdFile objects to combine
    :returns: the combined NrgFile object
    """
    # TODO: More sanity checks for agreeing parameters in the list
    if len(nrglist) == 1:
        return nrglist[0]
    result = nrglist.pop(0)
    for nrg in nrglist:
        # When times overlap, give following Nrgdata preference
        nrg.timearray = nrg.timearray.tolist()
        result.timearray = result.timearray.tolist()
        try:
            while result.timearray[-1] >= nrg.timearray[0]:
                del result.timearray[-1]
        except IndexError:
            raise IndexError("The nrg files completely overlap")
        result.dataarray = result.dataarray[:len(result.timearray), ...]
        result.timearray = np.concatenate((result.timearray, nrg.timearray))
        result.dataarray = np.concatenate((result.dataarray, nrg.dataarray))
    result.endtime = nrglist[-1].endtime
    del nrglist
    return result
