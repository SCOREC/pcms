#!/usr/bin/env python3
""" Command line tool for the python slice diagnostic"""
import argparse
import os

import pydiag.data.slices as sli
import pydiag.utils.comm as comm

xrange = (None, None)
yrange = (None, None)
zrange = (None, None)

parser = argparse.ArgumentParser()
parser.add_argument("time", help="Time window for averaging in the form start:end. "
                                 "Use f for the first and l for the last time in the data file.")
parser.add_argument("fileextension", nargs='+', help="List of file extensions of the runs"
                                                     " to be processed.")
parser.add_argument("--xavg", help="Average over full x direction", action="store_true")
parser.add_argument("--yavg", help="Average over full x direction", action="store_true")
parser.add_argument("--zavg", help="Average over full x direction", action="store_true")
parser.add_argument("--runpath", default=".", help="Path of the run")
parser.add_argument("--specname", help="Name of the species to process")
parser.add_argument("--var", default="phi",
                    help="Name of the variable, any of phi, apar, bpar, dens, tpar, tperp, qpar, "
                         "qperp, upar")

args = parser.parse_args()
os.chdir(args.runpath)

fileextension = args.fileextension[0]

if args.time.split(':')[0] == 'f':
    starttime = -1
elif args.time.split(':')[0] == 'l':
    starttime = -2
else:
    starttime = int(float(args.time.split(':')[0]))
if args.time.split(':')[1] == 'f':
    endtime = -1
elif args.time.split(':')[1] == 'l':
    endtime = -2
else:
    endtime = int(float(args.time.split(':')[1]))

common = comm.CommonData(fileextension, starttime, endtime)
diagspace = comm.DiagSpace(common.spatialgrid, x_fourier=common.x_local, y_fourier=common.y_local,
                           z_fourier=False, xrange=xrange, yrange=yrange, zrange=zrange, xavg=args.xavg,
                           yavg=args.yavg, zavg=args.zavg)
outputslice = sli.MomFieldSlice(common, args.var, args.specname, diagspace)

outputslice.generate_timeseries()
outputslice.write_to_numpy("slice_{}_{}".format(args.specname, fileextension))
outputslice.write_to_ascii("slice_{}_{}.dat".format(args.specname, fileextension))
