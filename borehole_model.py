# -*- coding: utf-8 -*-
"""
Created on Mon Jul 01 12:53:00 2019

@author: D. Loukrezis

Borehole model, retrieved from:
https://www.sfu.ca/~ssurjano/borehole.html
    
The Borehole function models water flow through a borehole. Its simplicity and 
quick evaluation makes it a commonly used function for testing a wide variety 
of methods in computer experiments. The response is water flow rate, in m3/yr.     
"""

import numpy as np

def borehole(x):
    """Analytical model of the water flow through a borehole"""
    rw = x[0] # radius of borehole (m)
    r  = x[1] # radius of influence (m)
    Tu = x[2] # transmissivity of upper aquifer (m^2/yr)
    Hu = x[3] # potentiometric head of upper aquifer (m)
    Tl = x[4] # transmissivity of lower aquifer (m^2/yr)
    Hl = x[5] # potentiometric head of lower aquifer (m)
    L  = x[6] # length of borehole (m)
    Kw = x[7] # hydraulic conductivity of borehole (m/yr)
    #
    frac_nom = 2 * np.pi * Tu * (Hu - Hl)
    frac_help = 2*L*Tu/(np.log(r/rw)*rw*rw*Kw)
    frac_denom = np.log(r/rw) * (1 + frac_help + Tu/Tl)
    #
    response = frac_nom / frac_denom
    return response