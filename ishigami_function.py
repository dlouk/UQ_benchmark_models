# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:58:19 2019

@author: loukrezis

Ishigami function
"""

import numpy as np

# function to be approximated
def ishigami(yvec):
    """Ishigami function"""
    a = 7
    b = 0.1
    y1 = yvec[0]
    y2 = yvec[1]
    y3 = yvec[2]
    term1 = np.sin(y1)
    term2 = a * np.sin(y2)**2
    term3 = b * y3**4 * np.sin(y1)
    return term1 + term2 + term3