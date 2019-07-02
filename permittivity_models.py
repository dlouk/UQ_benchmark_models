# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:00:30 2018

@author: D. Loukrezis

Permmitivity models for dielectric materials.
Uncertain parameters are typically modelled as uniform random variables.
The distribution ranges are given as: 
[nominal value - \alpha*(nominal value), nominal value + \alpha*(nominal value)].
Since relatively large uncertainties arise in these models, the value of \alpha 
is usually in the range [10%, 30%].
"""

import numpy as np

def debye1(freq, eps8, epss, tau, sigmas=0):
    """eps8: permittivity at high frequencies
    epss: permittivity in the static case
    tau: relaxation time"""
    eps0 = 8.8541878176 * 1e-12 # F/m
    omega = 2*np.pi*freq
    deps = epss - eps8
    epsr = eps8 + deps/(1.0 + omega*tau*1j) + sigmas/(omega*eps0*1j)
    return epsr

def cole1(freq, eps8, epss, tau, alpha, sigmas=0):
    eps0 = 8.8541878176 * 1e-12 # F/m
    omega = 2*np.pi*freq
    deps = epss - eps8
    epsr = eps8 + deps/(1.0 + (omega*tau*1j)**(1-alpha)) + \
           sigmas/(omega*eps0*1j)
    return epsr

def debye4(freq, eps8, epss1, epss2, epss3, epss4, tau1, tau2, tau3, tau4,
           sigmas=0):
    eps0 = 8.8541878176 * 1e-12 # F/m
    omega = 2*np.pi*freq
    term1 = (epss1 - eps8) / (1.0 + omega*tau1*1j)
    term2 = (epss2 - eps8) / (1.0 + omega*tau2*1j)
    term3 = (epss3 - eps8) / (1.0 + omega*tau3*1j)
    term4 = (epss4 - eps8) / (1.0 + omega*tau4*1j)
    epsr = eps8 + term1 + term2 + term3 + term4 + sigmas/(omega*eps0*1j)
    return epsr

def cole4(freq, eps8, epss1, epss2, epss3, epss4, tau1, tau2, tau3, tau4,
          alpha1, alpha2, alpha3, alpha4, sigmas=0):
    eps0 = 8.8541878176 * 1e-12 # F/m
    omega = 2*np.pi*freq
    term1 = (epss1 - eps8) / (1.0 + (omega*tau1*1j)**(1.0-alpha1))
    term2 = (epss2 - eps8) / (1.0 + (omega*tau2*1j)**(1.0-alpha2))
    term3 = (epss3 - eps8) / (1.0 + (omega*tau3*1j)**(1.0-alpha3))
    term4 = (epss4 - eps8) / (1.0 + (omega*tau4*1j)**(1.0-alpha4))
    epsr = eps8 + term1 + term2 + term3 + term4 + sigmas/(omega*eps0*1j)
    return epsr


if __name__ == "__main__":

	### Debye Test
	
    # Data from paper: "Modeling Human Tissues Using Fourth Order Debye Model
	# in Convolution-Based Three-Dimensional Finite-Difference Time-Domain",
	# Samah Mustafa, Amin M. Abbosh and Phong Thanh Nguyen,  
	# IEEE Transactions on Antennas and Propagation, 2014.
  
	# brain - white matter
    eps8_nom = 1.0
    epss1_nom = 33.009 + eps8_nom
    epss2_nom = 5.902 + eps8_nom
    epss3_nom = 45.645 + eps8_nom
    epss4_nom = 1.028*1e3 + eps8_nom
    tau1_nom = 8.261*1e-12
    tau2_nom = 127.4*1e-12
    tau3_nom = 1.795*1e-9
    tau4_nom = 21.134*1e-6
    sigmas_nom = 268.6*1e-3
    #
    epsr_d = debye4(1.0e9, eps8_nom, epss1_nom, epss2_nom, epss3_nom, epss4_nom,
                  tau1_nom, tau2_nom, tau3_nom, tau4_nom, sigmas_nom)
    print np.abs(epsr_d)
    
	
	### Cole-Cole Test
	
    # Data from paper: "The dielectric properties of biological tissues: III. 
	# Parametric models for the dielectric spectrum of tissues", 
	# S Gabriel, R W Lau and C Gabriel
	# Phys. Med. Biol., 1996.
	
    # brain - white matter
    eps8_nom = 4.0
    epss1_nom = 32.0 + eps8_nom
    epss2_nom = 100.0 + eps8_nom
    epss3_nom = 4.0*1e4 + eps8_nom
    epss4_nom = 3.5*1e7 + eps8_nom
    tau1_nom = 7.96*1e-12
    tau2_nom = 7.96*1e-9
    tau3_nom = 53.05*1e-6
    tau4_nom = 7.958*1e-3
    alpha1_nom = 0.1
    alpha2_nom = 0.1
    alpha3_nom = 0.3
    alpha4_nom = 0.02
    sigmas_nom = 0.02
    #
    epsr_c = cole4(1.0e9, eps8_nom, epss1_nom, epss2_nom, epss3_nom, epss4_nom,
                  tau1_nom, tau2_nom, tau3_nom, tau4_nom, alpha1_nom, alpha2_nom,
                  alpha3_nom, alpha4_nom,sigmas_nom)
    print np.abs(epsr_c)
