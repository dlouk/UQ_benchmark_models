# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:40:19 2017

@author: D. Loukrezis

Waveguide model with dielectric filling.
	- TE_10 excitation.
    - No dispersion.
    - Up to 7 input random variables.

Output: Reflection coefficient S_11.
"""

import numpy as np


def dispfree_single_freq(freq=10., width=30., height=3., fill_l=7., offset=5.,
                         epsr=2.0, muer=2.4, res='abs'):
    """Input parameters:
    1. frequency,
    2. width,
    3. height,
    4. filling length,
    5. vacuum offset,
    6. (constant w.r.t. frequency) permittivity,
    7. (constant w.r.t. frequency) permeability,

    Output: S11 parameter.
    res = 'dB' or 'abs' or 'complex' or 'im' or 're' """

    # scalings
    freq = freq * 1e9      # scale to GHz
    width = width * 1e-3   # scale to mm
    height = height * 1e-3 # scale to mm
    fill_l = fill_l * 1e-3 # scale to mm
    offset = offset * 1e-3 # scale to mm

    # constants
    c0 = 299792458 # speed of light in vacuum
    mue0 = 4 * np.pi * 1e-7 # vacuum permeability
    eps0 = 1.0 / (mue0 * c0**2.0) # vacuum permitivity

    # frequency
    w = 2*np.pi*freq # omega


    # wavenumbers
    k0 = w/c0 # in vacuum
    k = k0*np.sqrt(muer*epsr) # in material

    # propagation constants
    ky = np.pi / width # in y-direction
    kz = np.sqrt(k0**2 - ky**2) # in z-direction for vacuum
    kz1 = np.sqrt(k**2 - ky**2) # in z-direction for the debye-material

    # wave impedance for waveguide (TE-wave)
    Z0 = mue0*w/kz # for vacuum
    Z1 = mue0*muer*w/kz1 # for the filling

    # coefficient used to determine the reflection and transmission
    # coefficients from the equation system
    A = (-Z0/Z1 + 1)*np.exp(-1j*(offset + fill_l)*kz1)
    B = (Z0/Z1 + 1)*np.exp(1j*(offset + fill_l)*kz1)

    # reflection coefficient at the surface between the 2 materials
    r2 = (2*Z1*np.exp(-1j*offset*kz)) / ((Z1-Z0)*(np.exp(1j*kz1*offset) - \
                    np.exp(1j*(2*fill_l + offset)*kz1)*((Z1+Z0)/(Z1-Z0))**2))
    # transmission coefficient at the surface between the 2 materials
    t2 = -r2*B/A
    # reflection coefficient at the input port
    r1 = 0.5*np.exp(-1j*offset*kz)*(t2*np.exp(-1j*offset*kz1)*(-Z0/Z1 + 1) + \
                  r2*np.exp(1j*offset*kz1)*(Z0/Z1+1))

    # calculation of the S-parameters using the definition
    if res == 'complex':
        S11complex = r1
        return S11complex
    elif res == 'abs':
        S11abs = np.abs(r1)
        return S11abs
    elif res == 'dB':
        S11dB = 20*np.log10(np.abs(r1))
        return S11dB
    elif res == 'im':
        S11imag = r1.imag
        return S11imag
    elif res == 're':
        S11real = r1.real
        return S11real
    else:
        print "Not acceptable 'res' input."
        return


def dispfree_broadband(fmin=8., fmax=12., samples=1001, width=30.,
                              height=3., fill_l=10., offset=5., epsr=2.5,
                              muer=1., res='abs'):
    """Input parameters:
    1. frequency,
    2. width,
    3. height,
    4. filling length,
    5. vacuum offset,
    6. (constant w.r.t. frequency) permittivity,
    7. (constant w.r.t. frequency) permeability,

    Output: S11 parameter.
    res = 'dB' or 'abs' or 'complex' or 'im' or 're' """

    # get vector of operating frequencies
    freq = np.linspace(fmin, fmax, samples)

    # scalings
    freq = freq * 1e9      # scale to GHz
    width = width * 1e-3   # scale to mm
    height = height * 1e-3 # scale to mm
    fill_l = fill_l * 1e-3 # scale to mm
    offset = offset * 1e-3 # scale to mm

    # constants
    c0 = 299792458 # speed of light in vacuum
    mue0 = 4 * np.pi * 1e-7 # vacuum permeability
    eps0 = 1.0 / (mue0 * c0**2.0) # vacuum permitivity

    # frequency
    w = 2*np.pi*freq # omega


    # wavenumbers
    k0 = w/c0 # in vacuum
    k = k0*np.sqrt(muer*epsr) # in material

    # propagation constants
    ky = np.pi / width # in y-direction
    kz = np.sqrt(k0**2 - ky**2) # in z-direction for vacuum
    kz1 = np.sqrt(k**2 - ky**2) # in z-direction for the debye-material

    # wave impedance for waveguide (TE-wave)
    Z0 = mue0*w/kz # for vacuum
    Z1 = mue0*muer*w/kz1 # for the filling

    # coefficient used to determine the reflection and transmission
    # coefficients from the equation system
    A = (-Z0/Z1 + 1)*np.exp(-1j*(offset + fill_l)*kz1)
    B = (Z0/Z1 + 1)*np.exp(1j*(offset + fill_l)*kz1)

    # reflection coefficient at the surface between the 2 materials
    r2 = (2*Z1*np.exp(-1j*offset*kz)) / ((Z1-Z0)*(np.exp(1j*kz1*offset) - \
                    np.exp(1j*(2*fill_l + offset)*kz1)*((Z1+Z0)/(Z1-Z0))**2))
    # transmission coefficient at the surface between the 2 materials
    t2 = -r2*B/A
    # reflection coefficient at the input port
    r1 = 0.5*np.exp(-1j*offset*kz)*(t2*np.exp(-1j*offset*kz1)*(-Z0/Z1 + 1) + \
                  r2*np.exp(1j*offset*kz1)*(Z0/Z1+1))

    # calculation of the S-parameters using the definition
    if res == 'complex':
        S11complex = r1
        return S11complex
    elif res == 'abs':
        S11abs = np.abs(r1)
        return S11abs
    elif res == 'dB':
        S11dB = 20*np.log10(np.abs(r1))
        return S11dB
    elif res == 'im':
        S11imag = r1.imag
        return S11imag
    elif res == 're':
        S11real = r1.real
        return S11real
    else:
        print "Not acceptable 'res' input."
        return