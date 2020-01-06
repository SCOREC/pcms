""" This module contains a number of (semi-)analytic predictions from neoclassical theory"""

import numpy as np
import pydiag.utils.averages as av
import scipy.integrate as intg
import scipy.special as sps


def calc_potwidth(prd, ispec=1):
    """ Calculate an estimate for the potato width

    :param prd: ProfileData object of a run
    :param ispec: Index of the species we want, starting from 1
    :returns: potato width for the first species (should be main ions)
    """
    mspec = prd.cm.pars["mass{}".format(ispec)]
    potwidth = (4*prd.q[0]**2 * prd.cm.pnt.rhostar**2 * prd.cm.pnt.major_R * prd.cm.pnt.minor_r**2 *
                2*mspec*prd.T0s[0])**(1.0/3)
    return potwidth


def chang_hinton(prd, species):
    """ Calculate the Chang-Hinton estimate for the neoclassical heat flux

    :param prd: ProfileData object of a run
    :param species: Index of species to use
    :returns: Radial profile (numpy array) with the C-H prediction
    """
    # TODO: Implement the more general case for multiple ion species (Z_eff != 1)
    Qcharr = []
    epsilon = np.array(prd.xs * prd.cm.pnt.minor_r / prd.cm.pnt.major_R)
    mspec = 1  # since we only support one ion species so far
    prspec = prd.prspec[species]
    # alpha = Zeff - 1
    Z = 1
    f1 = 1 + 1.5 * epsilon ** 2
    f2 = (1 - epsilon ** 2) ** 0.5
    try:
        coll = prd.cm.pnt.coll
    except AttributeError:
        coll = 0
    # for alpha == 0: mustar == nustar
    for tb in range(len(prspec.dataarray["Ts"])):
        # k1 = 0 # relevant only for Zeff != 1
        mustar = prspec.nustar[tb]
        T = prspec.dataarray["Ts"][tb]
        n = prspec.dataarray["ns"][tb]
        omt = prspec.dataarray["omts"][tb]
        k2 = ((0.66 + 1.88 * epsilon ** 0.5 - 1.54 * epsilon) /
              (1 + 1.03 * mustar ** 0.5 + 0.31 * mustar) * f1 +
              1.6 * 0.74 * mustar * epsilon /
              (1 + 0.74 * mustar * epsilon ** 1.5) * 0.5 * (f1 - f2))
        qchval = (16.0 / 3 / np.sqrt(np.pi) * n ** 2 *
                  prd.q ** 2 / epsilon ** 1.5 * (1 - epsilon ** 2) *
                  Z ** 4 * coll * np.sqrt(2 * T * mspec) * k2 * omt)
        Qcharr.append(qchval)
    return np.array(Qcharr)


def fb_param_hinton_haz(prd):
    """ Calculate the force balance parameter k based on Hinton, Hazeltine

    :param prd: Profile data object
    :returns: numpy array of k for each radial position
    """
    avnustar = av.mytrapz(prd.prspec[0].nustar, prd.prspec[0].timearray)
    eps = np.array(prd.xs * prd.cm.pnt.minor_r / prd.cm.pnt.major_R)
    hhpred1 = (1.17 - 0.35 * avnustar ** (1 / 2)) / (1 + 0.7 * avnustar ** (1 / 2) * eps ** 3)
    hhpred2 = (hhpred1 - 2.1 * avnustar ** 2 * eps ** 3) / (1 + avnustar ** 2 * eps ** 3)
    return hhpred2


def fb_param_hirsch_sig(prd):
    """ Calculate the force balance parameter k based on Hirschman, Sigmar, 1981 and Satake 2010
    :param prd: Profile data object
    :returns: numpy array of k for each radial position
    """
    eps = np.array(prd.xs * prd.cm.pnt.minor_r / prd.cm.pnt.major_R)
    ft = 1.96*np.sqrt(eps)
    fc = 1-ft
    nustar = av.mytrapz(prd.prspec[0].nustar, prd.prspec[0].timearray)
    k = np.zeros(prd.cm.pnt.nx0)
    for i in range(prd.cm.pnt.nx0):
        K11, eK11 = intg.quad(unintK, 0, 5, args=(1, 1, nustar[i], eps[i]))
        K12, eK12 = intg.quad(unintK, 0, 5, args=(1, 2, nustar[i], eps[i]))
        K13, eK13 = intg.quad(unintK, 0, 5, args=(1, 3, nustar[i], eps[i]))
        mu1 = K11
        mu2 = K12 - 2.5*K11
        mu3 = K13 - 5*K12 + 25.0/4*K11
        k[i] = -np.sqrt(2)*mu2/(np.sqrt(2)*mu1 + ft[i]/fc[i]*(mu1*mu3-mu2**2))
    return np.array(k)


def unintK(x, i, j, nustar, eps):
    """ The matrix elements K_ij before x integration """
    nustar = np.array(nustar)
    return 8/3/np.sqrt(np.pi)*x**4*np.exp(-x**2)*x**(2*(i+j-2))*nutottau(x, nustar, eps)


def nutottau(x, nustar, eps):
    """ nu_tot(x)*tau_ii"""
    c = 3*np.sqrt(np.pi)/4
    return c*nu_d_nuhat(x)/(1+nustar*nu_d_nuhat(x)*c/x)/(1+15*np.pi**(3/2)/32*nustar*eps*nu_t_omT(x)/x)


def nu_d_nuhat(x):
    return (sps.erf(x)-chandr(x))/x**3


def nu_t_omT(x):
    return (sps.erf(x)-3*chandr(x))/x**3 + 8*chandr(x)/x


def chandr(x):
    """ The Chandrasehkar function"""
    return (sps.erf(x) - x*2/np.sqrt(np.pi)*np.exp(-x**2))/2/x**2
