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
            mpl.rc('font', size=15)
        else:
            self.cmap_unidirect = mpl.cm.hot
            mpl.rc('font', size=20)
        # for quantities with critical value (e.g. 0 for phi)
        self.cmap_bidirect = mpl.cm.bwr
        self.color_list = plt.cm.Dark2(np.linspace(0, 1.0, 9))  # Plot line colors
        # Set some basic plot properties
        # Those are set globally everytime a Plotting object is instantiated!
        mpl.rc('legend', fontsize="small")
        mpl.rc('lines', linewidth=2, markeredgewidth=1.5, markersize=10)
        # This sets the figure frames to transparent both on screen and in calls to savefig
        mpl.rc('figure', facecolor=(1, 1, 1, 0))
        mpl.rc('figure', frameon=False)
        mpl.rc('savefig', facecolor=(1, 1, 1, 0))
        mpl.rc('savefig', frameon=False)
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
                       'Qturb': r'$Q_{\mathrm{turb}}/Q_{gB}$',
                       'Qnc': r'$Q_{\mathrm{NC}}/Q_{gB}$',
                       'Piturb': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$',
                       'Pinc': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$',
                       "jbs": r'$j_\mathrm{BS}$', "phi": r"$\phi$",
                       "apar": r"$A_\parallel$", "dens": r"$n$",
                       "tpar": r"$T_\parallel$", "tperp": r"$T_\perp$",
                       "qpar": r"$q_\parallel + 1.5p_0 u_\parallel$",
                       "qperp": r"$q_\perp + p_0 u_\parallel$", "upar": r"$u_\parallel$",
                       "bpar": r"$B_\parallel$", "n00": r"$N_{00}$",
                       "n20": r"$N_{20}$", "n02": r"$N_{02}$"}

    @staticmethod
    def _minlogaxis(arr):
        """ Calculate lower limit for plotting on a log axis

        lower Integer * 10^minimumexponent

        :param arr: Array to base the limit on
        """
        minim = arr[arr != 0].min()  # ignore 0
        min_exp = np.floor(np.log10(minim))
        min_mant = np.floor(minim*10**(-min_exp))
        floormin = min_mant * 10**min_exp
        return floormin

    @staticmethod
    def _maxlogaxis(arr):
        """ Calculate upper limit for plotting on a log axis

         higher Integer * 10^maximumexponent

        :param arr: Array to base the limit on
        """
        maxim = arr[arr != 0].max()
        max_exp = np.floor(np.log10(maxim))
        max_mant = np.ceil(maxim*10**(-max_exp))
        ceilmax = max_mant * 10**max_exp
        return ceilmax

    @staticmethod
    def gyrobohm_SI(common, quantity):
        """ Convert gyroBohm unit to SI unit

        :param common: CommonData object
        :param quantity: Name of the quantity to calculate gB for
        """
        pnt = common.pnt
        elementary_charge = 1.60217662e-19
        temp_unit = 1e3*elementary_charge  # temp is in keV
        dens_unit = 1e19  # m**(-3)
        mass_unit = 1.6726219e-27  # proton mass
        # Magnetic field and reference length are already in T and m
        nref = pnt.nref*dens_unit
        Tref = pnt.Tref*temp_unit
        mref = pnt.mref*mass_unit
        cref = np.sqrt(Tref/mref)
        rhoref = mref*cref/elementary_charge/pnt.Bref
        Gamma_gB = nref*cref*rhoref ** 2/pnt.Lref**2
        gBunits = {"Ts": Tref, "ns": nref, "omts": 1, "omns": 1, "Gammanc": Gamma_gB/1e19,
                   "Gammaturb": Gamma_gB/1e19, "Qturb": Gamma_gB*Tref/1e3, "Qnc": Gamma_gB*Tref/1e3,
                   "Pinc": Gamma_gB*cref*mref, "Piturb": Gamma_gB*cref*mref,
                   "jbs": pnt.Bref*nref*cref*rhoref/pnt.Lref}
        return gBunits[quantity]

    @staticmethod
    def round_to_n(x, n):
        # Helper function to round to n significant digits
        return np.round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))