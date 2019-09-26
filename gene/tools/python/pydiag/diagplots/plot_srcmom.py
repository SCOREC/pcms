""" plot_srcmom.py: Plotting for the binary files of the source moment diagnostic"""


import matplotlib.pyplot as plt
from pydiag.diagplots.baseplot import Plotting
from pydiag.data.srcmom_data import SrcmomSeries
from pydiag.utils.averages import mytrapz


class PlotSrcmomAverage(Plotting):
    """Plot time averages of the source profiles

    :param srcmomseries: List of SrcmomData objects we want to plot
    """

    def __init__(self, srcmomseries):
        super().__init__()

        # Lists for the plotted structures
        self.leglist = []
        self.figlist = None
        self.axlist = None
        # Data structures used in the different plots
        self.src = srcmomseries
        for src in self.src:
            src.generate_timeseries()

    def createfigs(self):
        numfig = len(self.src[0].srcmoms)
        self.figlist = [plt.figure() for _ in range(numfig)]
        self.axlist = [fig.add_subplot(111) for fig in self.figlist]
        for srcav in self.src:
            self.plot_srcmoment(self.axlist, srcav)
        for ax in self.axlist:
            ax.legend()

    def plot_srcmoment(self, axlist, srcav):
        for isrc, srcdat in enumerate(srcav.srcmoms):
            axlist[isrc].plot(srcav.cm.spatialgrid.x_a, mytrapz(srcdat, srcav.timearray),
                              label=srcav.momtype + " " + srcav.srctype)
            axlist[isrc].set_xlabel(r'$x/a$')
            axlist[isrc].set_title(srcav.cm.specnames[isrc])


class PlotSrcmomxt(Plotting):
    """ Plot x-t colormaps of the source moments"""

    def __init__(self, srcmomaverage):
        super().__init__()
