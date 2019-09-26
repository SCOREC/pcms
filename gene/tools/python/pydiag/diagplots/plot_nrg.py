import numpy as np
import pydiag.utils.averages as av
from pydiag.diagplots.baseplot import Plotting, plt
from pydiag.utils import errors as err


class PlotNrgdata(Plotting):
    """ PlotNrgdata: Class to generate plots from nrgdata objects

    :param nrgdata: List of NrgData objects to plot
    """
    def __init__(self, nrgdata):
        super().__init__()
        for nrgd in nrgdata:
            nrgd.generate_timeseries()
            if len(nrgd.dataarray) == 0:
                print("No nrg data present in object ", nrgd.cm.fileextension)
                # Names of the columns in LaTeX readable form
        self.colnames = [r"|n_1|^2", r"|u_{1,\parallel}|^2", r"|T_{1,\parallel}|^2",
                         r"|T_{1,\perp}|^2", r"\Gamma_{es}", r"\Gamma_{em}", r"Q_{es}",
                         r"Q_{em}", r"\Pi_{es}", r"\Pi_{em}"]
        self.nrgdata = nrgdata
        self.n_col1 = 4
        self.n_col2 = 6
        if self.n_col1+self.n_col2 != nrgdata[0].cm.pnt.nrgcols:
            raise RuntimeError("The number of nrg columns in the Nrgdata class does not agree "
                               "with the current implementation in the plotting routine")

    def print_mean_err(self):
        """ Calculate and print the mean and error estimate for each nrg column """
        for nrgd in self.nrgdata:
            mean_nrg, err_nrg, ctimes_nrg = self._calc_mean_and_error(nrgd)
            print("Mean and errors for ", nrgd.objectname)
            for i_sp in range(nrgd.cm.pnt.n_spec):
                print("Species: ", i_sp, nrgd.cm.specnames[i_sp])
                print("============================")
                for col in range(nrgd.cm.pnt.nrgcols):
                    print(self.colnames[col]+": ", mean_nrg[i_sp, col], "+-", err_nrg[i_sp, col])
                print("corr times: ", ctimes_nrg[i_sp, :])
                print("============================")

    @staticmethod
    def _calc_mean_and_error(nrgd):
        mean_nrg = av.mytrapz(nrgd.dataarray, nrgd.timearray)
        err_nrg, ctimes_nrg = err.windowerr(nrgd.dataarray, nrgd.timearray)
        return mean_nrg, err_nrg, ctimes_nrg

    def plot_ttrace(self, poutput=False):
        """ Plot nrg file time traces with means and errors

        :param poutput: Generate pdf files of the plots as well
        """
        is_nonlinear = self._is_nrglist_nonlinear()
        fig1 = plt.figure()
        fig2 = plt.figure()
        max_spec = max([nrgd.cm.pnt.n_spec for nrgd in self.nrgdata])
        for nrgd in self.nrgdata:
            for i_sp in range(max_spec):
                ax1 = fig1.add_subplot(1, max_spec, i_sp+1)
                ax1.set_xlabel(r'$t\: [L_\mathrm{ref}/c_\mathrm{ref}]$')
                ax1.set_title(nrgd.cm.specnames[i_sp])
                ax2 = fig2.add_subplot(1, max_spec, i_sp+1)
                ax2.set_xlabel(r'$t\: [L_\mathrm{ref}/c_\mathrm{ref}]$')
                ax2.set_title(nrgd.cm.specnames[i_sp])
                if is_nonlinear:
                    self.plot_ttrace_nonlinear(ax1, ax2, nrgd, i_sp)
                else:
                    self.plot_ttrace_linear(ax1, ax2, nrgd, i_sp)
        ax1.legend()
        ax2.legend()
        fig1.tight_layout()
        fig2.tight_layout()
        if poutput == 1:
            fig1.savefig("nrg_ttrace1{}.pdf".format(self.nrgdata[0].cm.fileextension))
            fig2.savefig("nrg_ttrace2{}.pdf".format(self.nrgdata[0].cm.fileextension))

    def _is_nrglist_nonlinear(self):
        if np.all([nrgd.cm.nonlinear for nrgd in self.nrgdata]):
            return True
        elif np.all([not nrgd.cm.nonlinear for nrgd in self.nrgdata]):
            return False
        else:
            raise ValueError("Inconsistent combination of linear and nonlinear runs in nrgdata")

    def plot_ttrace_linear(self, ax1, ax2, nrgd, i_sp):
        ax1.set_yscale('log')
        ax2.set_yscale('log')
        self._find_rescalings(nrgd, i_sp)
        for col in range(self.n_col1):
            texlabel = r"$" + self.colnames[col] + r"$"
            base_line, = ax1.plot(nrgd.timearray, nrgd.dataarray[:, i_sp, col],
                                  label=texlabel)
        for col in range(self.n_col1, self.n_col1 + self.n_col2):
            texlabel = r"$" + self.colnames[col] + r"$"
            ax2.plot(nrgd.timearray, nrgd.dataarray[:, i_sp, col], label=texlabel)
        ax1.set_ylim(bottom=1e-3)
        ax2.set_ylim(bottom=1e-3)
        self._estimate_growthrate(nrgd, i_sp)

    @staticmethod
    def _find_rescalings(nrgd, i_sp):
        """ Find rescalings in linear runs"""
        col = 3  # Tperp should be a good quantity for this
        trace = nrgd.dataarray[:,  i_sp, col]
        rel = (trace/np.roll(trace, -1))[:-1]  # Find a large jump down between neighbouring points
        # Ignore the last element of rel, since it is trace[0]/trace[-1]
        return np.where(np.logical_and(rel > 1e12, trace[:-1] > 0.05))[0]

    def _estimate_growthrate(self, nrgd, i_sp):
        resc = self._find_rescalings(nrgd, i_sp)
        regrcolumns = [2]
        if len(resc) >= 2:
            start = resc[-2]+1    # Take time between the last two rescalings
            end = resc[-1] - 2    # Shift index by one to avoid including the rescaling
        elif len(resc) == 1:
            start = resc[0]+1     # Take time trace after rescaling
            end = None
        else:
            start = 0             # Take full time trace
            end = None
        for icol in regrcolumns:
            trace = nrgd.dataarray[start:end, i_sp, icol]/nrgd.dataarray[int(start), i_sp, icol]
            time = nrgd.timearray[start:end]-nrgd.timearray[start]
            growthrate = (np.log(trace[-1])-np.log(trace[0]))/(time[-1]-time[0])
            print("Growthrate estimate for {} {}: {}".format(nrgd.cm.specnames[i_sp],
                                                             self.colnames[icol], growthrate))

    def plot_ttrace_nonlinear(self, ax1, ax2, nrgd, i_sp):
        mean_nrg, err_nrg, ctimes_nrg = self._calc_mean_and_error(nrgd)

        # Plot first group: density, temperature, parallel flow etc
        for col in range(self.n_col1):
            texlabel = r"$" + self.colnames[col] + r"$"
            base_line, = ax1.plot(nrgd.timearray, nrgd.dataarray[:, i_sp, col],
                                  label=texlabel)
            if err_nrg[i_sp, col] != 0.0:
                ax1.axhline(y=mean_nrg[i_sp, col], linestyle='--', color=base_line.get_color())
                ax1.axhspan(mean_nrg[i_sp, col] - err_nrg[i_sp, col],
                            mean_nrg[i_sp, col] + err_nrg[i_sp, col], alpha=0.3,
                            facecolor=base_line.get_color())
        # Plot second group: Fluxes
        for col in range(self.n_col1, self.n_col1 + self.n_col2):
            texlabel = r"$" + self.colnames[col] + r"$"
            base_line, = ax2.plot(nrgd.timearray, nrgd.dataarray[:, i_sp, col],
                                  label=texlabel)
            ax2.axhline(y=mean_nrg[i_sp, col], linestyle='--', color=base_line.get_color())
            ax2.axhspan(mean_nrg[i_sp, col] - err_nrg[i_sp, col],
                        mean_nrg[i_sp, col] + err_nrg[i_sp, col], alpha=0.3,
                        facecolor=base_line.get_color())
