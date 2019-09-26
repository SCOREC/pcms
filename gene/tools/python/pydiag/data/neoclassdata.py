"""
Module handling neoclass files

"""
import csv
import numpy as np

from pydiag.data.base_file import TimeSeries


class NeoclassFile(TimeSeries):
    """NeoclassFile: Class to read a neoclass file from GENE

    :param filename: Typically neoclass.dat
    :param common: CommonData object of the run
    """
    def __init__(self, filename, common):
        super().__init__(filename, common)
        self.dataarray = []
        self.n_col = 4  # Particle, heat, momentum flux and bootstrap current
        self.n_spec = common.pnt.n_spec

    def generate_timeseries(self):
        """ Fill the NeoclassFile object with data """
        n_col = self.n_col
        try:
            with open(self.objectname) as ncfile:
                csvnrg = csv.reader(ncfile, delimiter=' ', skipinitialspace=True)
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
            raise IOError("neoclass file does not exist or has"
                          " wrong number of columns: {}".format(self.objectname))
        self.check_times()
        pos = self.calc_positions()
        self.dataarray = np.array(self.dataarray).astype(float, copy=False)
        # Reduce the nrgcols to only the required time frame
        self.dataarray = self.dataarray[pos, ...]
        self.timearray = np.array(self.timearray)[pos]
