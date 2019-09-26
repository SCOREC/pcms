""" Module containing fourier transform routines for GENE time series"""

import numpy as np


def kx_to_x(var_kx, nx0, axis=-3):
    """ Perform inverse FFT on kx spectral direction of variable

    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    :param var_kx: Variable in kx space
    :param nx0: Number of kx grid points
    :param axis: Which axis of var_kx is the x direction, by default the third last one
    :returns: variable in real x space
    """
    var_x = np.fft.ifft(nx0*var_kx, axis=axis)
    var_x = np.real_if_close(var_x, tol=1e5)
    return var_x


def ky_to_y(var_ky, nky0, axis=-2):
    """ Perform inverse FFT on ky spectral direction of variable

    The GENE data only include the non-negative ky components, so we need to use the real
    valued FFT routines
    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    :param var_ky: Variable in kx space
    :param nky0: Number of ky grid points
    :param axis: Which axis of var_kx is the x direction, by default the third last one
    :returns: variable in real y space
    """
    # The GENE data only include the non-negative ky components, so we need to use the real
    # valued FFT routines
    var_y = np.fft.irfft(2*nky0*var_ky, n=2*nky0, axis=axis)
    var_y = np.real_if_close(var_y, tol=1e5)
    return var_y


def kz_to_z(var_kz, nz0, axis=-1):
    """ Perform inverse FFT on kz spectral direction of variable

    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    This does not apply to the z direction though, since there is no z-spectral
    version of GENE.
    :param var_kz: Variable in kz space
    :param nz0: Number of kz grid points
    :param axis: Which axis of var_kz is the z direction, by default the third last one
    :returns: variable in real z space
    """
    var_z = np.fft.ifft(var_kz, axis=axis)
    var_z = np.real_if_close(var_z, tol=1e5)
    return var_z


def x_to_kx(var_x, nx0, axis=-3):
    """ Perform FFT on x direction of variable

    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    :param var_x: Variable in real x space
    :param nx0: Number of x grid points
    :param axis: Which axis of var_kx is the x direction
    :returns: variable in fourier kx space
    """
    var_kx = np.fft.fft(var_x, axis=axis)
    return var_kx


def y_to_ky(var_y, nky0, axis=-2):
    """ Perform FFT on y direction of variable

    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    :param var_y: Variable in real y space
    :param nky0: Number of y grid points
    :param axis: Which axis of var_ky is the y direction
    :returns: variable in fourier ky space
    """
    var_ky = np.fft.rfft(var_y, axis=axis)
    return var_ky


def z_to_kz(var_z, nz0, axis=-1):
    """ Perform FFT on z direction of variable

    Note: The fft and ifft in python and GENE/IDL have the factor 1/N switched!
    :param var_z: Variable in real z space
    :param nz0: Number of z grid points
    :param axis: Which axis of var_kz is the z direction
    :returns: variable in fourier kz space
    """
    var_kz = np.fft.fft(var_z, axis=axis)
    return var_kz
