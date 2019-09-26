#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from bisect import bisect_left, bisect_right
from copy import deepcopy
from dataclasses import dataclass


class Loader(object):
    """ This does the nasty job. We need for a given timestep to be able to load
        any kind of data. The idea is the following:
            1- the user sets folder and runs
            2- GUI determines which files and variables are available
            3- GUI fetches all times from all the files
            4- user selects a time window and a series of diagnostics
            5- GUI takes the subset of files that are needed
            5- GUI takes the requested subinterval in time
            6- loop over all times

            -- I win.
            """

    def __init__(self, diags, run, data, t_start=None, t_end=None):
        self.t_start = t_start
        self.t_end = t_end
        """ we need to call each diagnostic, get the file it needs """
        self.specnames = run.specnames
        self.needed_vars = []
        for diag in diags:
            self.needed_vars.append(diag.set_options(run, None))

        """ we need to decide which files need to be loaded by joining needed_vars 
            We simply join the whole but keep the original list"""
        global_vars = deepcopy(self.needed_vars[0])
        for file_key in global_vars.keys():  # $lis(self.needed_vars):
            for quant_key in global_vars[file_key].keys():
                for i_diag_needs in list(self.needed_vars):
                    global_vars[file_key][quant_key] = global_vars[file_key][quant_key] or \
                                                       i_diag_needs[file_key][quant_key]

        """ now we set needed_field to true or false depending if any diag needs it """
        self.needed_files = {}
        #        TODO what if data is not available
        """ loop over available fields"""
        for k in data.avail_vars.keys():
            """ loop over requested fields """
            if k in global_vars.keys():
                if any(global_vars[k].values()):
                    self.needed_files[k] = True
                else:
                    self.needed_files[k] = False
            else:
                self.needed_files[k] = False
        return

    def set_interval(self, data, t_start, t_end):

        class Processable_Time(object):
            def __init__(self, time):
                self.time = time
                pass

        """ Here we intersect all times with the one we want  """
        all_times = np.array([])

        """ loop over files that are necessary """
        for file in self.needed_files:

            if self.needed_files[file]:
                if len(all_times) == 0:
                    all_times = getattr(getattr(data.avail_times, file), 'times')
                else:
                    all_times = np.intersect1d(all_times,
                                               getattr(getattr(data.avail_times, file), 'times'))

        """ now reduce the interval """
        i_s = bisect_left(all_times, t_start)
        i_e = bisect_right(all_times, t_end)
        if i_s == len(all_times):
            i_s -= 1
        all_times = np.array(all_times[i_s:i_e])
        self.t_start = all_times[0]
        self.t_end = all_times[-1]

        """ and now fill the step list"""
        self.processable_times = all_times
        self.processable = [Processable_Time(t) for t in all_times]
        for file in self.needed_files:
            if self.needed_files[file]:
                tmp_t = getattr(getattr(data.avail_times, file), 'times')
                ind = np.arange(tmp_t.shape[0])[np.in1d(tmp_t, all_times)]
                my_steps = getattr(getattr(data.avail_times, file), 'steps')[ind]
                my_exts = [getattr(getattr(data.avail_times, file), 'files')[i] for i in ind]
                [setattr(self.processable[x], file, self.__TimeStep__(my_steps[x], my_exts[x][:]))
                 for x in range(all_times.size)]

        """ NOTE this does a synch between all diags. meaning if there is diag1 that wants field
        and diag2 that wants field and mom, diag1 will be executed only on the times of diag2
        This is a loss of data. To fix this we need to keep track of whether a given diag wants 
        to use multiple
        files at once or not."""

    @dataclass
    class __TimeStep__:
        step: int
        file: str
