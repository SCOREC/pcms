#!/usr/bin/env python3
""" Command line tool for the python diagnostics"""

import argparse
import sys
import os

import pydiag.data.fieldlib
import pydiag.data.srcmom_data
import pydiag.data.nrgdata as nrg
import pydiag.diagplots.plot_field as pfd
import pydiag.diagplots.plot_fsamom as pfsa
import pydiag.diagplots.plotpdf as ppdf
import pydiag.utils.averages
import pydiag.utils.comm as comm
from pydiag.data import profile_data, fsamom_data
from pydiag.diagplots import plot_profiles, plot_nrg
from pydiag.diagplots import plot_srcmom

# This set controls what is plotted.
PLOTTHIS = [
    'plotT',
    'plotQar',
    'plotQ',
#    'plotch',
    # 'plotP',
    'plotj',
    'plotG',
    'plotturb',
    'plotnu',
    ]
PLOTTHIS = set(PLOTTHIS)


parser = argparse.ArgumentParser()
parser.add_argument("time", help="Time window for averaging in the form start:end. "
                    "Use f for the first and l for the last time in the data file.")
parser.add_argument("fileextension", nargs='+', help="List of file extensions of the runs"
                                                     " to be processed.")
parser.add_argument("--runpath", default=".", help="Path of the run")
parser.add_argument("--continuation", "-c", action='store_true',
                    help="Try to combine the given files into one averaged profile"
                    "This requires matching parameter sets")
parser.add_argument("--makepdf", "-p", action='store_true', help="Generate pdf files of the plots")
parser.add_argument("--noturb", "-t", action='store_true',
                    help="Exclude turbulent transport in plots")
parser.add_argument("--timetrace", action='store_true',
                    help="Generate text files with time traces(currently only T at 3 fixed x/a)")
parser.add_argument("--profxt", "-x", action='store_true',
                    help="Create x-t contour plots instead of time averaged profiles")
parser.add_argument("--averageprofiles", action='store_true', help="Generate a time averaged profile"
                    "suitable for GENE input")
parser.add_argument("--field", "-f", action='store_true', help="Plot average and x-t contour of phi"
                                                               "instead of anything else")
parser.add_argument("--nrg", "-n", action='store_true', help="Perform mean and error calculation"
                    "for nrg file")
parser.add_argument("--probabilitydens", "-d", action='store_true', help="Generate a plot of the "
                                                                    "probability density function")
parser.add_argument("--fsa", action='store_true', help="Plot time averaged plots of the FSA moment"
                                                       "diagnostic")
parser.add_argument("--srcmom", action="store_true", help="Plot time averaged source moments")

args = parser.parse_args()
os.chdir(args.runpath)

fileextension = args.fileextension

if args.noturb:
    PLOTTHIS.discard('plotturb')

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

inputcomm = []
inputdata = []
inputfsa = []
inputfield = []
inputsrcmom = []

for fext in fileextension:
    common = comm.CommonData(fext, starttime, endtime)
    inputcomm.append(common)

if args.nrg:
    for common in inputcomm:
        filename = 'nrg{}'.format(common.fileextension)
        nrgdat = nrg.NrgFile(filename, common)
        nrgdat.readnrg()
        inputdata.append(nrgdat)
    if args.continuation:
        inputdata = [nrg.gluenrgdata(inputdata)]
    plotn = plot_nrg.PlotNrgdata(inputdata)
    plotn.print_mean_err()
    plotn.plot_ttrace(poutput=args.makepdf)
    plotn.show()

elif args.field:
    for common in inputcomm:
        fielddat = pydiag.data.fieldlib.FieldAverage(common)
        fielddat.generate_timeseries()
        inputfield.append(fielddat)
    if args.continuation:
        inputfield = [pydiag.data.fieldlib.gluefieldaverage(inputfield)]
    zon = pfd.Zonal(inputcomm, inputfield)
    zon.createfigs(withext=True, poutput=args.makepdf)
    plotf = pfd.PlotField(inputcomm, inputfield)
    plotf.plot_phi_xt(poutput=args.makepdf)
    plotf.show()
elif args.fsa:
    for common in inputcomm:
        try:
            fsadat = fsamom_data.FSAmomData(common)
            fsadat.getfsadata()
        except IOError:
            fsadat = fsamom_data.FSAmomData_legacy(common)
            fsadat.getfsadata()
        inputfsa.append(fsadat)
    plotfsa = pfsa.PlotFSA(inputfsa)
    plotfsa.createfigs(poutput=args.makepdf)
    plotfsa.show()
elif args.srcmom:
    for common in inputcomm:
        srcmom_f0 = pydiag.data.srcmom_data.SrcmomSeries(common, "energy", "f0")
        srcmom_f0.generate_timeseries()
        srcmom_ckh = pydiag.data.srcmom_data.SrcmomSeries(common, "energy", "ck_heat")
        srcmom_ckh.generate_timeseries()
        srcmom = [srcmom_f0, srcmom_ckh]
        inputsrcmom.append(srcmom)
    inputsrcmom = [item for sublist in inputsrcmom for item in sublist]

    plotsrc = plot_srcmom.PlotSrcmomAverage(inputsrcmom)
    plotsrc.createfigs()
    plotsrc.show()
else:
    for common in inputcomm:
        prdat = profile_data.ProfileData(common)
        fielddat = pydiag.data.fieldlib.FieldAverage(common)
        prdat.get_profiles()
        inputdata.append(prdat)
        if not (args.averageprofiles or args.timetrace or args.probabilitydens):
            fielddat.generate_timeseries()
            inputfield.append(fielddat)
    if args.continuation:
        inputdata = [profile_data.glueprofiledata(inputdata)]
        if not (args.averageprofiles or args.timetrace):
            inputfield = [pydiag.data.fieldlib.gluefieldaverage(inputfield)]
    if args.timetrace:
        for prdat in inputdata:       # Time trace generation
            prdat.generate_timetrace(0.45, "Qturb")
            prdat.generate_timetrace(0.6, "Qturb")
            prdat.generate_timetrace(0.75, "Qturb")

    if args.profxt:
        plot = plot_profiles.PlotProfxt(inputdata)
        plot.createfigs(args.makepdf)
        plot.show()
    elif args.probabilitydens:
        plotpd = ppdf.PlotProfilePdf(inputdata)
        plotpd.createfigs(poutput=args.makepdf)
        plotpd.show()
    else:
        plot = plot_profiles.PlotProfdata(inputdata, inputfield, PLOTTHIS)
        if args.averageprofiles:
            plot.output_avprof("average")
            sys.exit("Only average profile generated")
        plot.createfigs(args.makepdf)
        plot.show()
