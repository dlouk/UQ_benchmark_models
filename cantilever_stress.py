# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:59:41 2019

@author: loukrezis

Stress response on a uniform cantilever beam of width w and thickness t due to
the horizontal load Ph and the vertical load Pv.
"""

def cantilever_stress(yvec):
    """Cantilever beam stress function"""
    w = yvec[0]
    t = yvec[1]
    Pv = yvec[2]
    Ph = yvec[3]
    stress = 600*Pv + 600*Ph/(w*t*t)
    return stress