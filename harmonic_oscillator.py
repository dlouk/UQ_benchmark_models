# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 11:44:12 2019

@author: D. Loukrezis


Harmonic oscillator model.
    - Up to 6 random variables.

Parameter ranges:
x1: [1.0 - 3*0.05, 1.0 + 3*0.05]
x2: [1.0 - 3*0.01, 1.0 + 3*0.1]
x3: [0.1 - 3*0.01, 0.1 + 3*0.01]
x4: [0.5 - 3*0.05, 0.5 + 3*0.05]
x5: [0.45 - 3*0.075, 0.45 + 3*0.075]
x6: [1.0 - 3*0.2, 1.0 + 3*0.2]
"""

import numpy as np

def harm_osc(x):
    """Analytical model of harmonic oscillator"""
    w0 = np.sqrt( (x[1] + x[2]) / x[0] )
    response = 3*x[3] - np.abs( 2*x[4]/(x[0]*w0*w0) * np.sin(0.5*w0*x[5]))
    return response