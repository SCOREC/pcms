""" Base classes for file reading"""

import struct
from bisect import bisect_left, bisect_right
from os.path import getsize
import numpy as np
import utils.averages as av
import h5py
from abc import abstractmethod


class File(object):
    """ Base class to read files from GENE runs"""

    @abstractmethod
    def __init__(self, file=None, run_data=None, varname=None, prepath=None, spec=None,
                 extensions=None):
        pass

    @abstractmethod
    def redirect(self, extension):
        """ Call this routine to read from a new file"""
        pass

    @classmethod
    def __add_method__(self, name, idx):

        def _method(self, Time=None, Step=None):
            if Time:
                self.time = Time
            if not self.loaded_time[idx] == Time:
                self.loaded_time[idx] = Time
                if Step:
                    self.tind = Step.step
                    if self.extension != Step.file:
                        self.extension = Step.file
                        self._redirect(self.extension)
                if not self.fid:
                    self._redirect()
                setattr(self, name + "3d", self._readvar(idx))

            return getattr(self, name + "3d")

        return _method

    """ The following methods are general, for any file since is operation on time 
    array"""

    def FileName(self, extension=None):
        if self.run_data:
            return self.run_data.folder + '/' + (
                self.file + '_' + self.spec if self.spec else self.file) + (
                       extension if extension else self.extension)
        else:
            return None

    def set_times_and_inds(self):
        """ this is to be used for collecting all times from a series of extensions
            will return times, steps, and file extension.
            the whole time series is also set in the object"""
        self.reset_tinds()

        time_list = self.get_timearray(self.extensions[0])
        step_list = np.array(range(0, len(time_list), 1)).tolist()
        file_list = [self.extensions[0]]*len(time_list)

        for ext in self.extensions[1:]:
            t_loc = self.get_timearray(ext)
            time_list = np.concatenate((time_list, t_loc))
            step_list = np.concatenate((step_list, np.array(range(0, len(t_loc), 1))))
            file_list = file_list + [ext]*len(time_list)

        """ remove duplicates"""
        self.timeseries, inds = np.unique(time_list, return_index=True)
        time_list = np.array([time_list[i] for i in inds])
        step_list = np.array([step_list[i] for i in inds])
        file_list = [file_list[i] for i in inds]

        self.reset_tinds()

        return time_list, step_list, file_list

    def get_minmaxtime(self):
        if not self.timearray:
            self.get_timearray()
        return self.timearray[0], self.timearray[-1]

    def get_timearray(self, sparsefactor=1):
        raise NotImplementedError("Need to implement a method to read times")

    def _readvar(self, *args):
        raise NotImplementedError("Need to implement a method to read data")

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
        self.time = self._find_nearest_time(time)
        self.tind = self.timearray.index(self.time)

    def set_tind(self, tind):
        """ Set current timestep by index"""
        self.tind = tind
        self.time = self.timearray[tind]

    def get_var(self, name):
        """ """
        varidx = {v: k for k, v in self.varname.items()}
        if not self.loaded_time[varidx[name]] == self.time:
            self.loaded_time[varidx[name]] = self.time
            setattr(self, name + "3d", self._readvar(varidx[name]))
        return getattr(self, name + "3d")

    def reset_tinds(self):
        """ Reset time indices when reloading file"""
        self.loaded_time = [-1]*len(self.varname)
        self.tind = -1
        self.time = None
        self.timearray = None
        self.extension = self.extensions[0]
        self.filename = self.FileName()
        self._redirect(self.extension)

    def my_vars(self, to_print=None):
        if to_print:
            print(self.varname)
        return (self.varname)


class BinaryFile(File):
    """ Base class to read Fortran binary (unformatted) files from GENE runs"""

    def __init__(self, file=None, run_data=None, varname=None, prepath=None, spec=None,
                 nfields=None, extensions=None):
        """ I dont see the point of copying the run into the object,
        unless we want to create all objects at initialization and use them as list"""
        self.run_data = run_data
        self.file = file
        self.fid = None
        self.timearray = []
        self.tind = -1
        self.time = None
        self.nfields = nfields
        self.spec = spec
        self.extensions = extensions
        self.extension = run_data.fileextension
        self.filename = self.FileName(self.extension)
        self.varname = varname
        self.loaded_time = [-1]*len(self.varname)
        self._set_datatypes()
        self.te, self.tesize = self.time_entry()
        self._set_sizes()

        for idx in self.varname:
            setattr(self, varname[idx] + "3d", None)

        self.loaded_tind = [-1]*len(self.varname)

        for idx in varname:
            new_method = File.__add_method__(varname[idx], idx)
            setattr(File, varname[idx], new_method)

    def _redirect(self, extension=None):
        """ Call this routine to read from a new file"""
        if extension:
            self.extension = extension
            self.filename = self.FileName(self.extension)
        try:
            self.fid.close()
        except (AttributeError, OSError):
            pass
        self.fid = open(self.filename, 'rb')
        self.timearray = []

    #        self.get_timearray()
    #        self.reset_tinds()

    def get_timearray(self, extension=None):
        """Get time array from file """
        # am I pointing to the right file?
        if extension and self.extension != extension:
            self._redirect(extension)
            self.extension = extension
        # if not pointing to anything
        if not self.fid:
            self._redirect()

        self.timearray = []
        self.fid.seek(0)
        for i in range(int(getsize(self.filename)/(self.leapfld + self.tesize))):
            self.timearray.append(float(self.te.unpack(self.fid.read(self.tesize))[1]))
            self.fid.seek(self.leapfld, 1)

        return self.timearray

    def _set_datatypes(self):
        try:
            self.bigendian = self.run_data.pars['ENDIANNESS'] == 'BIG'
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
            self.realsize = 8 if self.run_data.pars['PRECISION'] == 'DOUBLE' else 4
        except KeyError:
            self.realsize = 8
        self.complexsize = 2*self.realsize
        self.entrysize = self.run_data.pnt.nx0*self.run_data.pnt.nky0*self.run_data.pnt.nz0*self\
            .complexsize
        # jump in bytes in field files
        self.leapfld = self.nfields*(self.entrysize + 2*self.intsize)

    def time_entry(self):
        """ Defines the struct for a time entry """
        # Fortran writes records as sizeof(entry) entry sizeof(entry)
        if self.bigendian:
            timeentry = struct.Struct('>idi')
        else:
            timeentry = struct.Struct('=idi')
        return timeentry, timeentry.size

    def offset(self, var):
        """Calculate offset in field file for a given self.time and variable"""
        if var in [i for i in range(self.nfields)]:
            return self.tesize + self.tind*(self.tesize + self.leapfld) + var*(
                    self.entrysize + 2*self.intsize) + self.intsize

    def _readvar(self, var):
        """ Return 3d field data at the time set in self.time"""
        self.fid.seek(self.offset(var))
        var3d = np.fromfile(self.fid,
                            count=self.run_data.pnt.nx0*self.run_data.pnt.nky0*self.run_data.pnt
                            .nz0,
                            dtype=self.npct)
        # Bring array into x, y. z order
        if self.run_data.x_local and not self.run_data.y_local:  # y-global has yx order
            var3d = var3d.reshape(self.run_data.pnt.nky0, self.run_data.pnt.nx0,
                                  self.run_data.pnt.nz0, order="F")
            var3d = np.swapaxes(var3d, 0, 1)
        else:
            var3d = var3d.reshape(self.run_data.pnt.nx0, self.run_data.pnt.nky0,
                                  self.run_data.pnt.nz0, order="F")
        return var3d


class H5File(File):
    """ Base class to read HDF5 files from GENE runs"""

    def __init__(self, file=None, run_data=None, varname=None, prepath=None, spec=None,
                 nfields=None, extensions=None):
        self.run_data = run_data
        self.file = file
        self.timearray = []
        self.tind = -1
        self.time = None
        self.nfields = nfields
        self.spec = spec
        self.fid = None
        self.extensions = extensions
        self.extension = run_data.fileextension
        self.filename = self.FileName(run_data.fileextension if run_data else None)
        self.varname = varname
        self.prepath = prepath

        for idx in self.varname:
            setattr(self, varname[idx] + "3d", None)

        self.loaded_tind = [-1]*len(self.varname)
        self.loaded_time = [-1]*len(self.varname)

        for idx in varname:
            new_method = File.__add_method__(self.varname[idx], idx)
            setattr(File, varname[idx], new_method)

    def _redirect(self, extension=None):
        """ Call this routine to read from a new file"""
        if extension:
            self.extension = extension
        self.filename = self.FileName(self.extension)
        try:
            self.h5file.close()
        except (AttributeError, OSError):
            pass
        self.fid = h5py.File(self.filename, 'r')
        self.timearray = []

    #        self.get_timearray()
    #        self.reset_tinds()

    def get_timearray(self, extension=None):
        # am I pointing to the right file?
        if extension and self.extension != extension:
            self._redirect(extension)
            self.extension = extension
        if not self.fid:
            self._redirect()
        """Get time array for field file """
        self.timearray = self.fid.get(self.prepath + "time").value.tolist()
        return self.timearray

    def _readvar(self, var):
        """ Return 3d field data at the time set in self.time"""
        out = self.fid.get(
            self.prepath + self.varname[var] + "/" + '{:010d}'.format(self.tind)).value
        if len(out.dtype) == 2:
            out = out['real'] + 1.0j*out['imaginary']
        out = np.swapaxes(out, 0, 2)
        if self.run_data.x_local and not self.run_data.y_local:  # y-global has yx order
            out = np.swapaxes(out, 0, 1)
        else:
            pass
        return out


class GENE_file(object):
    @classmethod
    def create(cls, file_type, run_data, spec=None, extensions=None):
        if run_data.is_h5:
            return H5File(file=file_type, run_data=run_data,
                          varname=cls._Variables(file_type, run_data),
                          prepath=cls._Prepath(file_type, run_data, spec), spec=spec,
                          nfields=cls._Nfields(file_type, run_data), extensions=extensions)
        elif run_data.is_adios:
            raise NotImplementedError('ADIOS not yet implemented')
        else:
            return BinaryFile(file=file_type, run_data=run_data,
                              varname=cls._Variables(file_type, run_data),
                              prepath=cls._Prepath(file_type, run_data, spec), spec=spec,
                              nfields=cls._Nfields(file_type, run_data), extensions=extensions)

    @classmethod
    def _Variables(cls, file_type, run_data):
        VARS_TO_FILE_MAP = {'field': {0: 'phi', 1: 'A_par', 2: 'B_par'}, 'mom': (
            {0: "n", 1: "u_par", 2: "T_par", 3: "T_per", 4: "Q_es", 5: "Q_em", 6: 'Gamma_es',
             7: 'Gamma_em'} if run_data.is3d else {0: "dens", 1: "T_par", 2: "T_perp", 3: "q_par",
                                                   4: "q_perp", 5: "u_par", 6: 'densI1',
                                                   7: 'TparI1', 8: 'TppI1'}), }

        if file_type not in VARS_TO_FILE_MAP:
            raise ValueError('Bad file type {}'.format(file_type))

        VARS_TO_FILE_MAP = cls._check_dict(VARS_TO_FILE_MAP, run_data)

        return VARS_TO_FILE_MAP[file_type]

    @classmethod
    def _Prepath(cls, file_type, run_data, spec):
        PREPATH_TO_FILE_MAP = {'field': '/field/', 'mom': '/mom_{}/'.format(spec)}

        if file_type not in PREPATH_TO_FILE_MAP:
            raise ValueError('Bad file type {}'.format(file_type))

        return PREPATH_TO_FILE_MAP[file_type]

    @classmethod
    def _Nfields(cls, file_type, run_data):
        NFIELDS_TO_FILE_MAP = {'field': int(run_data.pars['n_fields']),
            'mom': int(run_data.pars['n_moms'])
            #                TODOdoes this inherits the trap_passing splitting?
        }

        if file_type not in NFIELDS_TO_FILE_MAP:
            raise ValueError('Bad file type {}'.format(file_type))

        return NFIELDS_TO_FILE_MAP[file_type]

    @staticmethod
    def _check_dict(dic, run_data):
        """ check if we really have EM parts"""
        if not run_data.electromagnetic:
            dic['field'] = {0: 'phi'}
            if not run_data.is3d:
                dic['mom'] = dict((k, (dic['mom'])[k]) for k in range(6))
        elif not run_data.bpar:
            dic['field'] = {0: 'phi', 1: 'A_par'}
            dic['mom'] = dict((k, dic['mom'][k]) for k in range(6))
        #        TODO add trap-passing splitting
        return dic


class TimeSeries(object):
    """ Base class for for any time series

    This is used for formatted (ascii) GENE output as well as derived data
    from GENE diagnostics
    :param objectname: Either a file name (e.g. profile_Ions.dat) or a (freely choosable) name for
     a derived dataset
    :param common: CommonData object of the current run
    """

    def __init__(self, objectname, run_data):
        self.objectname = objectname
        self.cm = run_data
        #        self.starttime = common.starttime
        #        self.endtime = common.endtime
        #        self.fileobject = None  # A binary file object providing data to this one
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
        if len(self.timearray) != len(set(self.timearray)):
            raise RuntimeError("Error: {} contains 2 blocks with "
                               "identical timestamp".format(self.objectname))
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

    def calc_nearest(self, target):
        """ For a single time, find element closest to given input time"""
        timearray = self.timearray
        if target >= timearray[-1]:
            return len(timearray) - 1
        elif target <= timearray[0]:
            return 0
        else:
            return bisect_right(timearray, target)


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
