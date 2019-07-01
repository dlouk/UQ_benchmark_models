# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 12:49:18 2019

@author: D. Loukrezis

Steel column cost-reliability function, retrieved from 
https://www.sfu.ca/~ssurjano/steelcol.html
"""

import numpy as np


def steel_column(x):
    """Steel column reliability-cost function"""
    Fs = x[0] # yield stress (MPa)
    Pd = x[1] # dead weight load (N)
    P1 = x[2] # variable load 1 (N)
    P2 = x[3] # variable load 2 (N)
    B  = x[4] # flange breadth (mm)
    D  = x[5] # flange thickness (mm)
    H  = x[6] # profile height (mm)
    F0 = x[7] # initial deflection (mm)
    E  = x[8] # Young's modulus (MPa)
    L  = x[9] # column length
    #
    P = Pd + P1 + P2 # total load
    Eb = (np.pi*np.pi)*E*B*D*H*H/(2*L*L) # Euler buckling load
    #
    frac1 = 1./(2*B*D)
    frac2 = (F0*Eb)/(B*D*H*(Eb-P))
    #
    response = Fs - P*(frac1 + frac2)
    return response