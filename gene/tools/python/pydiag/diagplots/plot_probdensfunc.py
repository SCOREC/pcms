"""
Module for plotting of probability density functions
"""

import numpy as np
from matplotlib.colors import LogNorm
from pydiag.diagplots.baseplot import plt, mpl, Plotting


class PlotProfilePdf(Plotting):
    """ Plot probability density distributions

    :param profdata: ProfileData object list to use
    """

    def __init__(self, profdata):
        super().__init__()
        # Lists for the plotted structures
        self.figlist = []
        self.axlist = []
        self.leglist = []
        # Data structures used in the different plots
        self.prdlist = profdata

    def createfigs(self, poutput):
        """ Create the canvas of figures and axes objects

        :param poutput: Switch to produce output as pdf files
        """
        ntypes = 2  # Number of different plot forms
        nspecs = len(self.prdlist[0].prspec)
        numfigs = nspecs * ntypes
        self.figlist = [plt.figure() for _ in range(numfigs)]
        # One plot window for each species and quantity, one axis object in each
        self.axlist = [fig.add_subplot(111) for fig in self.figlist]
        for iprd, prd in enumerate(self.prdlist):
            toplot = self.axlist[0:nspecs]
            self.plot_pdf_atx((0.4, 0.75), toplot, prd)
            self.plot_pdf_cont(self.axlist[nspecs:2*nspecs], prd)
        self.axlist[0].legend()
        for fig in self.figlist:
            fig.tight_layout()
        if poutput:
            self.figlist[0].savefig("pdfat06{}.pdf".format(self.prdlist[0].fext))

    def plot_pdf_atx(self, xpos, axlist, prd, quantity="Qturb"):
        """ Plot the pdf of a quantity at a certain radial position

        :param xpos: Radial position in fraction of the minor radius
        :param axlist: Axes objects to plot into
        :param prd: ProfileData object to plot from
        :param quantity: string which of the profile quantities should be used.
        Can be any of ns, Ts, omn, omt, Gammaturb, Gammanc, Qturb, Qnc, Piturb, Pinc, jbs
        """

        xmin = (np.abs(prd.xs - xpos[0])).argmin()
        xmax = (np.abs(prd.xs - xpos[1])).argmin()
        for iprs, prs in enumerate(prd.prspec):
            ttrace = np.array(prs.profdata[quantity])[:, xmin:xmax].flatten()
            pdf, binedges = np.histogram(ttrace, bins=50, range=(0, 7), density=True)
            linest = "-."
            try:
                if prd.pnt.include_f0_contr:
                    labstring = "w/ NC"
                    linest = "-"
            except AttributeError:
                labstring = "w/o NC"
            axlist[iprs].plot(binedges[:-1] + binedges[1] / 2, pdf, linest, label=labstring)
            axlist[iprs].set_xlabel(self.titles[quantity])

    def plot_pdf_cont(self, axlist, prd, quantity="Qturb"):
        """ Colormap of the the pdf

        :param axlist: Axes objects to plot into
        :param prd: ProfileData object to plot from
        :param quantity: string which of the profile quantities should be used.
        :returns: List of plotted colormaps (for colorbar)
        """
        for iprs, prs in enumerate(prd.prspec):
            arrquant = np.array(prs.profdata[quantity])
            minim = np.amin(arrquant)
            maxim = np.amax(arrquant)
            pdf = []
            for ix in range(len(prd.xs)):
                xpdf, bins = np.histogram(arrquant[:, ix], range=(minim, maxim),
                                          bins=50, density=True)
                pdf.append(xpdf)
            axlist[iprs].pcolormesh(prd.xs, bins[:-1], np.array(pdf).T, cmap=self.cmap_unidirect,
                                    norm=LogNorm(vmin=0.1, vmax=1))
            axlist[iprs].set_ylabel(self.titles[quantity])
            axlist[iprs].set_xlabel(r"$x/a$")
