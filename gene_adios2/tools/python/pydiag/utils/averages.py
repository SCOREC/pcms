import numpy as np


def mytrapz(yvar, timefld):
    """ Trapezoid rule, adjusted for single timestep

    Operates on the first dimension of yvar (typically time)
    hence the name timefld for the integration variable samples
    :param yvar: Variable to integrate over
    :param timefld: The integration variable (time)
    :returns: The integration result (dim(yvar)-1)
    """
    timefld = np.array(timefld)
    yvar = np.array(yvar, dtype=float)
    if timefld.size == 1:
        return np.atleast_1d(yvar)[0]
    else:
        if yvar.shape[0] != len(timefld):
            raise ValueError("First dimension of yvar and timefld do not match")
        tmpt = timefld / (timefld[-1] - timefld[0])
        return np.trapz(yvar, x=tmpt, axis=0)


def z_av3d(var, geom):
    """ Perform the average in z direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object

    """
    if geom.cm.x_local and geom.cm.y_local:
        return np.average(var, weights=geom.jacobian, axis=-1)
    else:
        return np.average(var, weights=geom.jaco3d.T, axis=-1)  # geom has (z, x or y) arrays


def y_av3d(var, geom):
    """ Perform the average in y direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """
    if geom.cm.y_local:   # Add the negative half of the Fourier space for ky
        var[:, 0, :] *= 0.5
        return np.sum(2*var, axis=1)
    else:
        return np.average(var, weights=geom.jaco3d.T, axis=1)


def x_av3d(var, geom):
    """ Perform the average in x direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """
    if geom.cm.x_local:
        if geom.cm.y_local:
            return np.sum(var, axis=0)
        else:   # Mirror positive kx to negative kx
            var[0, :, :] *= 0.5
            return np.sum(2*var, axis=0)
    else:
        return np.average(var, weights=geom.jaco3d.T, axis=0)


def xz_av3d(var, geom):
    """ Perform the average in x and z direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """
    if geom.cm.x_local:
        if geom.cm.y_local:
            return np.average(np.sum(var, axis=0), weights=geom.jacobian, axis=-1)
        else:   # Mirror positive kx to negative kx
            var[0, :, :] *= 0.5
            return np.average(np.sum(2*var, axis=0), weights=geom.jacobian.T, axis=-1)
    else:
        return np.average(var, weights=geom.jaco3d.T, axis=(0, -1))


def yz_av3d(var, geom):
    """ Perform the average in y and z direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """

    if geom.cm.y_local:   # Add the negative half of the Fourier space for ky
        var[:, 0, :] *= 0.5
        # The transposition of the Jacobian only matters for x-global
        return np.average(np.sum(2*var, axis=1), weights=geom.jacobian.T, axis=-1)
    else:
        return np.average(var, weights=geom.jaco3d.T, axis=(1, -1))  # geom has (z, x or y) arrays


def xy_av3d(var, geom):
    """ Perform the average in x and y direction for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """
    if geom.cm.x_local:
        if geom.cm.y_local:
            var[:, 0, :] *= 0.5
            return np.sum(2*var, axis=(0, 1))
        else:   # Mirror positive kx to negative kx
            var[0, :, :] *= 0.5
            return np.average(np.sum(2*var, axis=0), weights=geom.jacobian.T, axis=1)
    else:
        if geom.cm.y_local:
            var[:, 0, :] *= 0.5
            return np.average(np.sum(2*var, axis=1), weights=geom.jacobian.T, axis=0)
        else:
            return np.average(var, weights=geom.jaco3d.T, axis=(0, 1))


def xyz_av3d(var, geom):
    """ Perform the average in all spatial directions for a 3d variable

    :param var: Variable to average over
    :param geom: GENE geomtry object
    """
    if geom.cm.x_local:
        if geom.cm.y_local:
            var[:, 0, :] *= 0.5
            return np.average(np.sum(2*var, axis=(0, 1)), weights=geom.jacobian)
        else:   # Mirror positive kx to negative kx
            var[0, :, :] *= 0.5
            temp = np.sum(2*var, axis=0)
            return np.average(temp, weights=geom.jacobian.T, axis=(0, -1))
    else:
        if geom.cm.y_local:
            var[:, 0, :] *= 0.5
            return np.average(np.sum(2*var, axis=1), weights=geom.jacobian.T, axis=(0, -1))
        else:
            return np.average(var, weights=geom.jaco3d.T)


def av3d_by_switch(xavg, yavg, zavg):
    """ Map the required averaging function by binary switches for the averages

    :returns: Reference to the appropriate function from this module
    """
    if not (xavg or yavg or zavg):
        return lambda var, geom: var  # Return an identity function
    if xavg and not (yavg or zavg):
        return x_av3d
    if xavg and yavg and not zavg:
        return xy_av3d
    if xavg and not yavg and zavg:
        return xz_av3d
    if xavg and yavg and zavg:
        return xyz_av3d
    if not (xavg or yavg) and zavg:
        return z_av3d
    if not xavg and yavg and zavg:
        return yz_av3d
    if not xavg and yavg and not zavg:
        return y_av3d
