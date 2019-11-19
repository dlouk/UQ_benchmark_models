# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:49:03 2019

@author: loukrezis

Meromorphic function / N-dimensional
"""

import numpy as np

def qoi_mero_Nd(yvec):
    """Meromorphic function"""
    NN = len(yvec)
    gvec_tilde = np.linspace(1, NN, NN)
    gvec_tilde = 1.0/gvec_tilde**2
    coeff = 1.0 / (2.0*np.linalg.norm(gvec_tilde, ord=1))
    gvec = gvec_tilde * coeff
    dotprod = np.dot(gvec, yvec)
    return 1.0/(1 + dotprod)