usage: pydiagc.py [-h] [--runpath RUNPATH] [--continuation] [--xavg] [--yavg]
                  [--zavg] [--xfourier] [--yfourier] [--zfourier]
                  [--species SPECIES [SPECIES ...]] [--kyind KYIND] [--save]
                  [--timetrace] [--averageprofiles] [--gkdb] [--noturb]
                  [--profxt | --zonal | --nrg | --probabilitydens | --fsa | --srcmom | --spectra | --perpspectra | --fluxspectra | --xycontours | --ballooning | --anisotropy]
                  time fileextension [fileextension ...]

positional arguments:
  time                  Time window for averaging in the form start:end. Use f
                        for the first and l for the last time in the data
                        file.
  fileextension         List of file extensions of the runs to be processed.

optional arguments:
  -h, --help            show this help message and exit
  --runpath RUNPATH     Path of the run
  --continuation, -c    Try to combine time traces into one objectThis
                        requires matching parameter sets
  --noturb, -t          Exclude turbulent transport in profile plots
  --profxt, -x          Create x-t contour plots of profiles
  --zonal, -z           Plot average and x-t contour of zonal phiinstead of
                        anything else
  --nrg, -n             Plot nrg file with average and uncertainty estimate
  --probabilitydens, -d
                        Generate a plot of the probability density function
                        for fluxes in profiles
  --fsa                 Plot time averaged plots of the FSA moment diagnostic
  --srcmom              Plot time averaged source moments
  --spectra             Plot amplitude spectra of mom and field
  --perpspectra         Plot 1d amplitude spectra of fields vs. perpendicular
                        wavenumber
  --fluxspectra         Plot flux spectra based on mom and field
  --xycontours          Contour plots in x and y based on mom and field
  --ballooning, -ball   Plot mode structure for all fields
  --anisotropy          Plot scale-dependent anisotropy of magnetic field
                        fluctuations

Diagnostic grid options:
  --xavg                Average over full x direction
  --yavg                Average over full y direction
  --zavg                Average over full z direction
  --xfourier            Operate in kx space
  --yfourier            Operate in ky space
  --zfourier            Operate in kz space
  --species SPECIES [SPECIES ...], -s SPECIES [SPECIES ...]
                        Names of species to process

Index values / ranges:
  --kyind KYIND         ky index to be considered in diagnostics (if
                        supported)

Output options:
  --save                Save the plots as pdfs or videos (for animations)
  --timetrace           Generate text files with time traces (currently only T
                        at 3 fixed x/a) from profiles
  --averageprofiles     Generate a time averaged profile suitable for GENE
                        input
  --gkdb                Save data for GKDB (github.com/gkdb/gkdb)
