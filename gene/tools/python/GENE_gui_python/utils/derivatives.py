#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:20:41 2018

@author: merlo
"""
import numpy as np
from scipy.sparse import spdiags


def compute(dx, F, order, accuracy=None):
    F = np.reshape(F, (len(F), 1))

    if order == 1:
        return derivativeO1(dx, F, accuracy if accuracy else 2)
    elif order == 2:
        return derivativeO2(dx, F, accuracy if accuracy else 3)
    else:
        raise NotImplementedError("Not implemented derivatives higher than second order")


def derivativeO1(dx, F, accuracy):
    m = len(F)
    e = np.ones((1, m))
    z = np.zeros((1, m))
    if accuracy == 2:
        D = spdiags(np.row_stack((-0.5*e, z, 0.5*e)), np.arange(-1, 2), m, m)
    elif accuracy == 4:
        D = spdiags(np.row_stack((1/12*e, -2/3*e, z, 2/3*e, -1/12*e)), np.arange(-2, 3), m, m)
    elif accuracy == 6:
        D = spdiags(np.row_stack((-1/60*e, 3/20*e, -3/4*e, z, 3/4*e, -3/20*e, 1/60*e)),
                    np.arange(-3, 4), m, m)
    elif accuracy == 8:
        D = spdiags(
            np.row_stack((1/280*e, -4/105*e, 1/5*e, -4/5*e, z, 4/5*e, -1/5*e, 4/105*e, -1/280*e)),
            np.arange(-4, 5), m, m)
    else:
        raise NotImplementedError("Not implemented accuracy other than 2, 4, 6 or 8")

    dF = D*F
    dF[0] = F[1] - F[0]
    dF[-1] = F[-1] - F[-2]

    return dF/dx


def derivativeO2(dx, F, accuracy):
    m = len(F)
    e = np.ones((m, 1))
    if accuracy != 3:
        raise NotImplementedError("Not implemented accuracy other than 3")

    D = spdiags(np.row_stack((e, -2*e, e)), np.arange(-1, 2), m, m)

    dF = D*F;
    dF[0] = dF[1];
    dF[-1] = dF[-2]

    return dF/dx
