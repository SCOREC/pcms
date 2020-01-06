""" Base module for plotting

Contains:
- Matplotib imports and settings

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from pkg_resources import parse_version

if parse_version(mpl.__version__) >= parse_version("1.5"):
    from cycler import cycler


class Plotting(object):
    """ Base class for plotting routines

    Contains the parameters used for the plots such as
    line colors, and font sizes
    Sets the line and point marker widths/sizes
    Sets the default color cycle (respecting an API change in matplotlib 1.5)
    """

    def __init__(self):
        # Matplotlib 2.x changes the default plot styles
        if parse_version(mpl.__version__) >= parse_version("2.0"):
            # Colormap for quantities >0 (e.g. heat flux)
            self.cmap_unidirect = mpl.cm.inferno
            self.font = FontProperties(size=12)
            self.xyfs = 15
            self.titlefs = 21
        else:
            self.cmap_unidirect = mpl.cm.hot
            self.font = FontProperties(size=16)
            self.xyfs = 20
            self.titlefs = 28
        # for quantities with critical value (e.g. 0 for phi)
        self.cmap_bidirect = mpl.cm.bwr
        self.color_list = plt.cm.Dark2(np.linspace(0, 1.0, 9))  # Plot line colors
        # Set some basic plot properties
        # Those are set globally everytime a Plotting object is instantiated!
        mpl.rc('lines', linewidth=2, markeredgewidth=1.5, markersize=10)
        mpl.rc('xtick', labelsize=self.xyfs)
        mpl.rc('ytick', labelsize=self.xyfs)
        # mpl 1.5 introduces a new, more flexible prop_cycle parameter, so different line styles
        # can be defined in the else case as well
        if parse_version(mpl.__version__) <= parse_version("1.4"):
            mpl.rc('axes', color_cycle=self.color_list)
        else:
            mpl.rc('axes', prop_cycle=cycler('color', self.color_list))
        # Dictionary that converts internal names to well-formatted LaTeX for the plots
        self.titles = {'omts': r'$a/L_T$', 'omns': r'$a/L_n$', 'ns': r'$n$', 'Ts': r'$T$',
                       'Gammaturb': r'$\Gamma_{\mathrm{turb}}/\Gamma_{gB}$',
                       'Gammanc': r'$\Gamma_{\mathrm{NC}}/\Gamma_{gB}$',
                       'Qturb': r'$Q_{\mathrm{turb}}/Q_{gB}$', 'Qnc': r'$Q_{\mathrm{NC}}/Q_{gB}$',
                       'Piturb': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$',
                       'Pinc': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$', "jbs": r'$j_\mathrm{BS}$'}
