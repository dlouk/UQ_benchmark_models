# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 11:57:10 2019

@author: D. Loukrezis

Meromorphic function 1/(G*Y), where
G = [g1, g2, ..., g16 = ][1, 0.5, 0.1, 0.05, 0.01, ..., 5*1e-8]
Y = [Y1, Y2, ..., Y16], Y_n ~ U[-1,1]

Taken from: "Adaptive Polynomial Approximation by Means of Random Discrete 
Least Squares", G. Migliorati, ENUMATH 2013.
"""

import numpy as np

def mero(yvec):
    """Meromorphic function"""
    gvec_tilde = np.array([1e0, 5*1e-1, 1e-1, 5*1e-2, 1e-2, 5*1e-3, 1e-3,
                           5*1e-4, 1e-4, 5*1e-5, 1e-5, 5*1e-6, 1e-6, 5*1e-7,
                           1e-7, 5*1e-8])
    coeff = 1.0 / (2.0*np.linalg.norm(gvec_tilde, ord=1))
    gvec = gvec_tilde * coeff
    dotprod = np.dot(gvec, yvec)
    return 1.0/(1 + dotprod)