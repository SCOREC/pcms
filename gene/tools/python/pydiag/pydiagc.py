#!/usr/bin/env python3
""" Command line tool for the python diagnostics"""

import argparse
import sys
import os
import matplotlib.pyplot as plt
from multiprocessing.dummy import Pool

import pydiag.data.datafiles as datafiles
import pydiag.data.srcmom_data
import pydiag.data.nrgdata as nrg
import pydiag.diagplots.plot_zonal as pzf
import pydiag.diagplots.plot_fsamom as pfsa
import pydiag.diagplots.plot_probdensfunc as ppdf
import pydiag.utils.averages
import pydiag.utils.comm as comm
import pydiag.utils.gkdb as gkdb
import pydiag.diagplots.plot_spectra as spectra
import pydiag.diagplots.plot_contour as contour
import pydiag.diagplots.plot_torus as torus
import pydiag.diagplots.plot_ball as ballooning
import pydiag.diagplots.plot_anisotropy as anisotropy
from pydiag.data import profile_data, fsamom_data
from pydiag.diagplots import plot_profiles, plot_nrg
from pydiag.diagplots import plot_srcmom
from pydiag.data.base_file import gluetimetraces

xrange = (None, None)
yrange = (None, None)
zrange = (None, None)

# This set controls what is plotted. TODO: Get rid of this control structure
PLOTTHIS = ['plotT', 'plotQar', 'plotQ',  # 'plotch',
            # 'plotP',
            'plotj', 'plotG', 'plotturb', 'plotnu', ]
PLOTTHIS = set(PLOTTHIS)

parser = argparse.ArgumentParser()
# Necessary parameters
parser.add_argument("time", help="Time window for averaging in the form start:end. "
                                 "Use f for the first and l for the last time in the data file.")
parser.add_argument("fileextension", nargs='+', help="List of file extensions of the runs"
                                                     " to be processed.")

# Options related to file handling
parser.add_argument("--runpath", default=".", help="Path of the run")
parser.add_argument("--continuation", "-c", action='store_true',
                    help="Try to combine time traces into one object"
                         "This requires matching parameter sets")

# Options for the diagnostic grid to use
gridopts = parser.add_argument_group("Diagnostic grid options")
gridopts.add_argument("--xavg", help="Average over full x direction", action="store_true")
gridopts.add_argument("--yavg", help="Average over full y direction", action="store_true")
gridopts.add_argument("--zavg", help="Average over full z direction", action="store_true")
gridopts.add_argument("--xfourier", help="Operate in kx space", action="store_true")
gridopts.add_argument("--yfourier", help="Operate in ky space", action="store_true")
gridopts.add_argument("--zfourier", help="Operate in kz space", action="store_true")
gridopts.add_argument("--species", "-s", help="Names of species to process", nargs='+')

# Auxiliary options to request particular index values (/ranges)
indexargs = parser.add_argument_group("Index values / ranges")
indexargs.add_argument("--kyind", type=int, nargs=1, #action="store_const",
                       default=None,
                       help="ky index to be considered in diagnostics (if supported)")


# Output options
outputargs = parser.add_argument_group("Output options")
outputargs.add_argument("--save", action='store_true',
                        help="Save the plots as pdfs or videos (for animations)")
outputargs.add_argument("--timetrace", action='store_true',
                        help="Generate text files with time traces "
                             "(currently only T at 3 fixed x/a) from profiles")
outputargs.add_argument("--averageprofiles", action='store_true',
                        help="Generate a time averaged profile suitable for GENE input")
outputargs.add_argument("--gkdb", action='store_true',
                        help="Save data for GKDB (github.com/gkdb/gkdb)")

# Different types of diagnostics
parser.add_argument("--noturb", "-t", action='store_true',
                    help="Exclude turbulent transport in profile plots")
diagopts = parser.add_mutually_exclusive_group()
diagopts.add_argument("--profxt", "-x", action='store_true',
                      help="Create x-t contour plots of profiles")

diagopts.add_argument("--zonal", "-z", action='store_true',
                      help="Plot average and x-t contour of zonal phi"
                           "instead of anything else")
diagopts.add_argument("--nrg", "-n", action='store_true', help="Plot nrg file with average and "
                                                               "uncertainty estimate")
diagopts.add_argument("--probabilitydens", "-d", action='store_true',
                      help="Generate a plot of the probability density function for fluxes "
                           "in profiles")
diagopts.add_argument("--fsa", action='store_true', help="Plot time averaged plots of the"
                                                         " FSA moment diagnostic")
diagopts.add_argument("--srcmom", action="store_true", help="Plot time averaged source moments")
diagopts.add_argument("--spectra", action="store_true",
                      help="Plot amplitude spectra of mom and field")
diagopts.add_argument("--perpspectra", action="store_true",
                      help="Plot 1d amplitude spectra of fields vs. perpendicular wavenumber")
diagopts.add_argument("--fluxspectra", action="store_true",
                      help="Plot flux spectra based on mom and field")
diagopts.add_argument("--xycontours", action="store_true",
                      help="Contour plots in x and y based on mom and field")
diagopts.add_argument("--ballooning", "-ball", action="store_true",
                      help="Plot mode structure for all fields")
diagopts.add_argument("--anisotropy", action='store_true', help="Plot scale-dependent anisotropy "
                      "of magnetic field fluctuations")
diagopts.add_argument("--toruscut", action='store_true', help="Plot contours on"
                                                              " a cut through the torus")

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
    starttime = float(args.time.split(':')[0])
if args.time.split(':')[1] == 'f':
    endtime = -1
elif args.time.split(':')[1] == 'l':
    endtime = -2
else:
    endtime = float(args.time.split(':')[1])

inputcomm = []

inputdata = []
inputdiagspace = []
inputfsa = []
inputfield = []
inputmom = []
inputsrcmom = []
inputnrg = []
inputprof = []

speclist = []
for fext in fileextension:
    common = comm.CommonData(fext, starttime, endtime)
    inputcomm.append(common)
    if not args.species:
        speclist += common.specnames
    else:
        speclist = args.species
    diagspace = comm.DiagSpace(common.spatialgrid, x_fourier=args.xfourier, y_fourier=args.yfourier,
                               z_fourier=args.zfourier, xrange=xrange, yrange=yrange, zrange=zrange,
                               xavg=args.xavg, yavg=args.yavg, zavg=args.zavg)
    inputdiagspace.append(diagspace)
    rundatafiles = datafiles.RunDataFiles(common)
    inputdata.append(rundatafiles)

# Remove all double entries from species list
speclist = list(set(speclist))

if args.nrg:
    for icomm, common in enumerate(inputcomm):
        filename = 'nrg{}'.format(common.fileextension)
        nrgdat = inputdata[icomm].get_fileobject("nrg")
        inputnrg.append(nrgdat)
    if args.continuation:
        inputnrg = [nrg.gluenrgdata(inputnrg)]
    plotn = plot_nrg.PlotNrgdata(inputnrg)
    plotn.print_mean_err()
    plotn.plot_ttrace(poutput=args.save)


elif args.zonal:
    for icomm, common in enumerate(inputcomm):
        fielddat = pzf.ZonalAverage(common, inputdata[icomm], xrange=xrange)
        inputfield.append(fielddat)
    if args.continuation:
        inputfield = [gluetimetraces(inputfield)]
    zon = pzf.PlotZonalAverage(inputcomm, inputfield)
    zon.createfigs(withext=True, poutput=args.save)
    plotf = pzf.PlotZonalxt(inputcomm, inputfield)
    plotf.plot_phi_xt(poutput=args.save)
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
    plotfsa.createfigs(poutput=args.save)
elif args.srcmom:
    for icomm, common in enumerate(inputcomm):
        srcmom_f0 = pydiag.data.srcmom_data.SrcmomSeries(common, inputdata[icomm], "energy", "f0")
        srcmom_ckh = pydiag.data.srcmom_data.SrcmomSeries(common, inputdata[icomm], "energy",
                                                          "ck_heat")
        srcmom = [srcmom_f0, srcmom_ckh]
        inputsrcmom.append(srcmom)
    inputsrcmom = [item for sublist in inputsrcmom for item in sublist]

    plotsrc = plot_srcmom.PlotSrcmomAverage(inputsrcmom)
    plotsrc.createfigs()

elif args.spectra:
    inputampspec = []
    for icomm, common in enumerate(inputcomm):
        with Pool(len(speclist)) as p:
            inputampspec = p.map(
                    lambda spec: spectra.AmplitudeSpectra(common, spec, inputdata[icomm]),
                    speclist)
    specplt = spectra.PlotAmplitudeSpectra(inputcomm, inputampspec)
    specplt.createfigs()

elif args.perpspectra:
    inputperpspec = []
    for icomm, common in enumerate(inputcomm):
        with Pool(len(args.species)) as p:
            inputperpspec = p.map(
                    lambda spec: spectra.PerpSpectra(common, spec, inputdata[icomm]),
                    args.species)
    specplt = spectra.PlotPerpSpectra(inputcomm, inputperpspec)
    specplt.createfigs()

elif args.fluxspectra:
    inputfluxspec = []
    for icomm, common in enumerate(inputcomm):
        with Pool(len(speclist)) as p:
            inputfluxspec = p.map(lambda spec: spectra.FluxSpectra(common, spec, inputdata[icomm],
                                                                   with_momentum=False),
                                  speclist)
    specplt = spectra.PlotFluxSpectra(inputcomm, inputfluxspec)
    specplt.createfigs()
    for s in inputfluxspec:
        s.print_total_fluxes()

elif args.xycontours:
    inputcontours = []
    for icomm, common in enumerate(inputcomm):
        with Pool(len(speclist)) as p:
            inputcontours = p.map(
                    lambda spec: contour.MomFieldContour(common, spec, inputdata[icomm],
                                                         zavg=args.zavg,
                                                         quantities=["phi"]), speclist)
    contplot = contour.PlotXYContour(inputcontours)
    contplot.createfigs(animate=True, save=args.save)

elif args.toruscut:
    inputtorus = []
    for icomm, common in enumerate(inputcomm):
        spec = speclist[0]
        if len(speclist) != 1:
            print("Using only first given species for torus cut: {}".format(spec))
        inputtorus.append(torus.TorusContourCut(common, spec, inputdata[icomm], quantity="phi"))
    torusplot = torus.PlotTorusCut(inputtorus)
    torusplot.createfigs(animate=True, save=args.save)

elif args.ballooning:
    for icomm, common in enumerate(inputcomm):
        fielddat = ballooning.ModeStructure(common, inputdata[icomm], kyind=args.kyind,normalize=True)
        inputfield.append(fielddat)
    if args.continuation:
        inputfield = [gluetimetraces(inputfield)]
    ballplot = ballooning.PlotModeStructure(inputcomm, inputfield)
    ballplot.createfigs()

elif args.gkdb:
    for icomm, common in enumerate(inputcomm):
        gkdb_obj = gkdb.GKDB(common, inputdata[icomm])
        inputfield.append(gkdb_obj)
    if args.continuation:
        inputfield = [gluetimetraces(inputfield)]
    gkdb_obj.Write_json('./')
elif args.anisotropy:
    for icomm, common in enumerate(inputcomm):
        fielddat = anisotropy.Anisotropy(common, inputdata[icomm])
        inputfield.append(fielddat)
#    if args.continuation:
#        inputfield = [gluetimetraces(inputfield)]
    aniplot = anisotropy.PlotAnisotropy(inputcomm, inputfield)
    aniplot.createfigs()

else:
    for icomm, common in enumerate(inputcomm):
        prdat = profile_data.ProfileData(common, inputdata[icomm])
        fielddat = pydiag.diagplots.plot_zonal.ZonalAverage(common, inputdata[icomm], xrange=xrange)
        prdat.get_profiles()
        inputprof.append(prdat)
        if not (args.averageprofiles or args.timetrace or args.probabilitydens or args.profxt):
            fielddat.generate_timeseries()
            inputfield.append(fielddat)
    if args.continuation:
        inputprof = [profile_data.glueprofiledata(inputprof)]
        if not (args.averageprofiles or args.timetrace):
            inputfield = [gluetimetraces(inputfield)]
    if args.timetrace:
        for prdat in inputprof:  # Time trace generation
            prdat.generate_timetrace(0.45, "Qturb")
            prdat.generate_timetrace(0.6, "Qturb")
            prdat.generate_timetrace(0.75, "Qturb")

    if args.profxt:
        plot = plot_profiles.PlotProfxt(inputprof)
        plot.createfigs(args.save)

    elif args.probabilitydens:
        plotpd = ppdf.PlotProfilePdf(inputprof)
        plotpd.createfigs(poutput=args.save)

    else:
        plot = plot_profiles.PlotProfdata(inputprof, inputfield, PLOTTHIS)
        plot_profiles.printaverageflux(inputprof, 0.4, 0.6)
        if args.averageprofiles:
            plot.output_avprof("average")
            sys.exit("Only average profile generated")
        plot.createfigs(args.save)

plt.show()
