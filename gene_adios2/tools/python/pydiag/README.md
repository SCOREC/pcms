usage: pydiagc.py [-h] [--runpath RUNPATH] [--continuation] [--makepdf]
                  [--noturb] [--timetrace] [--profxt] [--averageprofiles]
                  [--field] [--nrg] [--probabilitydens] [--fsa] [--srcmom]
                  time fileextension [fileextension ...]

positional arguments:
  time                  Time window for averaging in the form start:end. Use f
                        for the first and l for the last time in the data
                        file.
  fileextension         List of file extensions of the runs to be processed.

optional arguments:
  -h, --help            show this help message and exit
  --runpath RUNPATH     Path of the run
  --continuation, -c    Try to combine the given files into one averaged
                        profileThis requires matching parameter sets
  --makepdf, -p         Generate pdf files of the plots
  --noturb, -t          Exclude turbulent transport in plots
  --timetrace           Generate text files with time traces(currently only T
                        at 3 fixed x/a)
  --profxt, -x          Create x-t contour plots instead of time averaged
                        profiles
  --averageprofiles     Generate a time averaged profilesuitable for GENE
                        input
  --field, -f           Plot average and x-t contour of phiinstead of anything
                        else
  --nrg, -n             Perform mean and error calculationfor nrg file
  --probabilitydens, -d
                        Generate a plot of the probability density function
  --fsa                 Plot time averaged plots of the FSA momentdiagnostic
  --srcmom              Plot time averaged source moments
