""" Base classes for file reading"""

import struct
from bisect import bisect_left, bisect_right
import mmap

import numpy as np
import pydiag.utils.averages as av


class BinaryFile(object):
    """ Base class to read Fortran binary (unformatted) files from GENE runs"""

    def __init__(self, filename, common, usemmap=False):
        self.cm = common
        self.file = None
        self.timearray = []
        self.tind = 0
        self.time = None
        self.filename = filename
        if usemmap:
            self.readfunc = np.frombuffer
        else:
            self.readfunc = self._npfromfilewithoffset
        self.usemmap = usemmap
        self._set_datatypes()
        self.te, self.tesize = self.time_entry()
        self._set_sizes()
        self.redirect(filename)  # Initialises self.file

    @staticmethod
    def _npfromfilewithoffset(file, count, dtype, offset):
        file.seek(offset)
        return np.fromfile(file, count=count, dtype=dtype)

    def redirect(self, filename):
        """ Call this routine to read from a new file"""
        self.filename = filename
        try:
            self.file.close()
        except (AttributeError, OSError):
            pass
        file = open(filename, 'rb')
        if self.usemmap:
            self.file = mmap.mmap(file.fileno(), length=0, access=mmap.ACCESS_READ)
        else:
            self.file = file
        self.timearray = []
        self.get_timearray()
        self.reset_tinds()

    def _set_datatypes(self):
        try:
            self.bigendian = self.cm.pars['ENDIANNESS'] == 'BIG'
        except KeyError:
            self.bigendian = False
        if self.bigendian:
            self.nprt = (np.dtype(np.float64)).newbyteorder()
            self.npct = (np.dtype(np.complex128)).newbyteorder()
        else:
            self.nprt = np.dtype(np.float64)
            self.npct = np.dtype(np.complex128)

    def _set_sizes(self):
        self.intsize = 4
        try:
            self.realsize = 8 if self.cm.pars['PRECISION'] == 'DOUBLE' else 4
        except KeyError:
            self.realsize = 8
        self.complexsize = 2*self.realsize

    def time_entry(self):
        """ Defines the struct for a time entry """
        # Fortran writes records as sizeof(entry) entry sizeof(entry)
        if self.bigendian:
            timeentry = struct.Struct('>idi')
        else:
            timeentry = struct.Struct('=idi')
        return timeentry, timeentry.size

    def _find_nearest_time(self, time):
        pos = bisect_left(self.timearray, time)
        if pos == 0:
            return self.timearray[0]
        if pos == len(self.timearray):
            return self.timearray[-1]
        before = self.timearray[pos - 1]
        after = self.timearray[pos]
        if after - time < time - before:
            return after
        else:
            return before

    def set_approximate_time(self, time):
        ti = self._find_nearest_time(time)
        self.set_time(ti)

    def set_time(self, time):
        """ Set current timestep by value"""
        self.time = time
        self.tind = self.timearray.index(time)

    def set_tind(self, tind):
        """ Set current timestep by index"""
        self.tind = tind
        self.time = self.timearray[tind]

    def get_minmaxtime(self):
        if not self.timearray:
            self.get_timearray()
        return self.timearray[0], self.timearray[-1]

    def get_timearray(self):
        raise NotImplementedError("Need to implement a method to read times from binary file")

    def reset_tinds(self):
        """ Reset time index and previously read data """
        self.tind = 0

    def offset(self, *args):
        raise NotImplementedError("Need to implement a method to calculate data offset in file")

    def readvar(self, *args):
        raise NotImplementedError("Need to implement a method to read binary data")

    def get_var(self, varname):
        raise NotImplementedError(
            "Need to implement a method to provide variable access by name {}".format(varname))


class TimeSeries(object):
    """ Base class for for any time series

    This is used for formatted (ascii) GENE output as well as derived data
    from GENE diagnostics
    :param objectname: Either a file name (e.g. profile_Ions.dat) or a (freely choosable) name for
     a derived dataset
    :param common: CommonData object of the current run
    """

    def __init__(self, objectname, common):
        self.objectname = objectname
        self.cm = common
        self.starttime = common.starttime
        self.endtime = common.endtime
        self.fileobject = None  # A binary file object providing data to this one
        self.timearray = []
        self.dataarray = []
        self.timeaverage = None

    def generate_timeseries(self):
        raise NotImplementedError("Need to implement a method to aggregate data")

    def generate_timeaverage(self):
        if len(self.dataarray):
            self.timeaverage = av.mytrapz(self.dataarray, self.timearray)
        else:
            self.generate_timeseries()
            self.timeaverage = av.mytrapz(self.dataarray, self.timearray)

    def get_minmaxtime(self):
        """ Return first and last time in the object

        If self.fileobject is set, try to get it from there
        """
        if np.array(self.timearray).size == 0:
            try:
                self.timearray = self.fileobject.timearray
            except AttributeError:
                raise RuntimeError("Time array not available or empty")
        return self.timearray[0], self.timearray[-1]

    def check_times(self):
        """ Set the boundaries of the time window to read"""
        first_time, last_time = self.get_minmaxtime()
       # if len(self.timearray) != len(set(self.timearray)):
        #    raise RuntimeError("Error: {} contains 2 blocks with "
         #                      "identical timestamp".format(self.objectname))
        if self.starttime == -1 or (0 < self.starttime < first_time):
            print("Using first time present in {} for starttime".format(self.objectname))
            self.starttime = first_time
        if self.endtime == -1:
            print("Using first time present in {} for endtime".format(self.objectname))
            self.endtime = first_time
        if self.starttime == -2:
            print("Using last time present in {} for starttime".format(self.objectname))
            self.starttime = last_time
        if (self.endtime == -2) or (self.endtime > last_time):
            print("Using last time present in {} for endtime".format(self.objectname))
            self.endtime = last_time
        if (self.endtime < first_time) or (self.starttime > last_time):
            raise RuntimeError("Time window not contained in {}".format(self.objectname))
        print(("starttime={}, endtime={},"
               " first_time={}, last_time={}".format(self.starttime, self.endtime, first_time,
                                                     last_time)))

    def calc_positions(self):
        """ For a single time, find element closest to given input time"""
        timearray = self.timearray
        if self.starttime == self.endtime:
            pos = np.array([bisect_left(timearray, self.starttime)])
        else:
            startpos = bisect_left(timearray, self.starttime)
            endpos = bisect_right(timearray, self.endtime)
            pos = np.arange(startpos, endpos)
        return pos


def gluetimetraces(tserieslist):
    """ Function to concatenate a list of time series (continuation runs) into one single object

    for continuation runs
    :param tserieslist: List of TimeSeries objects to connect
    :returns: One connected TimeSeries object
    """
    # TODO: Consistency checks if runs match
    if len(tserieslist) == 1:
        return tserieslist[0]
    # Use first TimeSeries as basis for the construction of the glued object
    result = tserieslist.pop(0)
    # If the timefield has been converted to a np array before we need to undo it
    try:
        result.timearray = result.timearray.tolist()
    except AttributeError:
        pass
    for tseries in tserieslist:
        # When times overlap, give following TimeSeries preference
        while result.timearray[-1] > tseries.timearray[0]:
            del result.timearray[-1]
            del result.dataarray[-1]
        result.timearray = np.concatenate((result.timearray, tseries.timearray))
        result.dataarray = np.concatenate((result.dataarray, tseries.dataarray))
    result.endtime = tserieslist[-1].endtime
    del tserieslist
    return result
