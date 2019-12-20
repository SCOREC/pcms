""" errors.py: Collection of routines to calculate error estimates and other statistical tools

"""
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from pydiag.utils.averages import mytrapz


def autocorrtime(data, timefld):
    """ Calculate autocorrelation time of data array

    :param data: First dimension is considered the direction to correlate over
    :param timefld: The domain to correlate in
    :returns: 1/e time of the autocorrelation function in an array
    that has ndim(data)-1 dimensions
    """
    data = np.array(data, dtype=float)
    if len(timefld) <= 2:
        return np.zeros(data.shape[1:])
    if data.ndim == 1:
        return autocorrtime_1d(data, timefld)
    else:
        return autocorrtime_nd(data, timefld)


def autocorrtime_1d(data, timefld):
    """ Correlation time calculation for 1d arrays (only time)

    :param data: 1d array of the data to correlate
    :param timefld: 1d array of the correlation domain
    :returns: 1/e correlation time
    """
    meanvar = mytrapz(data, timefld)
    # The mean needs to be substracted to arrive at
    # the statistical definition of the autocorrelation
    data -= meanvar
    result = np.correlate(data, data, mode="full")
    result = result[int(result.size/2):]
    overe = np.squeeze(np.where(result > np.exp(-1)*result[0]))
    if overe.size <= 1:
        corrtime = 0
    else:
        cort_ind = np.array_split(overe, np.where(np.diff(overe) != 1)[0] + 1)[0][-1]
        corrtime = timefld[cort_ind] - timefld[0]
    return corrtime


def autocorrtime_nd(data, timefld):
    """ Correlation time calculation for >1d arrays

    :param data: array of the data to correlate
    :param timefld: array of the correlation domain
    :returns: 1/e correlation times in a numpy array
    """
    flat = data.reshape(data.shape[0], -1)
    pool = ThreadPool()
    corrtimes = np.array(pool.map(lambda ind: _autocorrtime_nd_helper(ind, flat, timefld, data),
                                  np.ndindex(data.shape[1:])))
    # This parallelises the following loop
    # corrtimes = np.zeros(data.shape[1:])
    # for ind in np.ndindex(data.shape[1:]):
    #     corrtimes[ind] = _autocorrtime_nd_helper(ind, flat, timefld, data)
    pool.close()
    pool.join()
    return corrtimes.reshape(data.shape[1:])


def _autocorrtime_nd_helper(ind, flat, timefld, data):
    var = np.squeeze(flat[:, np.ravel_multi_index(ind, data.shape[1:])])
    return autocorrtime_1d(var, timefld)


def windowerr(data, timefld, std_mean=True, n_win=0):
    """ Wrapper for 1d or multidimensional windowed error calculation

    :param data: array of the data to calculate the uncertainty
    :param timefld: array of the sampling times for data
    :param std_mean: switches to error of the mean with student-t correction for
    low n_win instead of standard deviation of the windowed sample
    :param n_win: Number of windows to group the data. 0 means automatic calculation based on the
    autocorrelation time (recommended default)
    :returns: The windowed error and the autocorrelation time.
    """
    data = np.array(data)
    if n_win <= 0:
        ctimes = autocorrtime(data, timefld)
    else:  # Assign a dummy array for the correlation times
        ctimes = np.zeros(data.shape[1:])
    if data.ndim == 1:
        return windowerr_1d(data, timefld, ctimes, std_mean, n_win), ctimes
    else:
        return windowerr_nd(data, timefld, ctimes, std_mean, n_win), ctimes


def windowerr_1d(data, timefld, corrt, std_mean=True, n_win=0):
    """ Standard deviation of the windowed mean values

    Tries to use window widths of 5 corrt, fallback are at least 2
    std_mean switches to error of the mean
    with student-t correction for low n_win

    """
    if len(timefld) <= 2:
        # print("Not enough times for error calculation")
        return 0
    total_dt = timefld[-1] - timefld[0]
    # If no window number is given, use correlation time input
    # start with window width = 5 tcorr; if fewer than 6, reduce width, min 2
    # same as in IDL diag, tools.pro, stat_avg_err,
    # except the minimal n_win is 10 there
    if n_win <= 0:
        if total_dt >= 30*corrt:
            if np.allclose(corrt, 0):
                return 0
            else:
                n_win = int(np.floor(total_dt/(5*corrt)))
        else:
            if total_dt >= 12*corrt:
                n_win = 6
            else:
                # print('Not enough windows')
                return 0
    win_inds = []
    wwidth = total_dt/n_win
    for iwin in range(n_win):
        win_inds.append(np.squeeze(np.where(
            (timefld - timefld[0] >= iwin*wwidth)&(timefld - timefld[0] < (iwin + 1)*wwidth))))
    means = np.empty(n_win)
    for iwin, window in enumerate(win_inds):
        if np.array(window).size < 1:
            print('windowerr_1d: windows with < 1 steps, aborting error computation')
            return 0
        means[iwin] = mytrapz(data[window], timefld[window])
    std_error = np.std(means, ddof=1)
    if std_mean:
        std_error *= np.sqrt(1./(n_win - 2))
    return std_error


def windowerr_nd(data, timefld, ctimes, std_mean=True, n_win=0):
    """ Standard error of the windowed mean value

    Function for multidimensional data which calls the 1d version repeatedly
    """
    std = np.zeros(data.shape[1:])
    flat = data.reshape(data.shape[0], -1)
    for ind in np.ndindex(data.shape[1:]):
        var = np.squeeze(flat[:, np.ravel_multi_index(ind, data.shape[1:])])
        std[ind] = windowerr_1d(var, timefld, ctimes[ind], std_mean, n_win)
    return std
