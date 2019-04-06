# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 11:59:35 2019

@author: Bram-Notebook
"""

import numpy as np
import matplotlib.pyplot as plt


def snells(theta_n,theta_z,n):
    """
    SNELLS determines the angle (with regard to the optical axis) after the
    encounter with the interface using the angle of the incoming ray, the
    angle of the normal to the optical axis and refractive indices.
    """
    n0,n1 = n
    theta_1 = np.arcsin(np.sin(theta_n)*n0/n1)+theta_z
    return theta_1

def propagate2sph(r_0,theta_0,d,D,R):
    """
    PROPAGATE2SPH determines the location (z_1, r_1) of a ray with initial
    height (r_0) and angle (theta_0) that hits a spherical interface of
    diameter D and radius R at a distance d. Also, the angle to the
    surface's normal (theta_n) and the angle of the normal to the optical
    axis are calculated.   
    """
    
    a = 1 + np.tan(theta_0)**2
    b = 2*r_0*np.tan(theta_0)-2*R-2*d
    c = r_0**2+2*d*R+d**2
    
    if R > 0:
        z_1 = np.min(np.roots([a,b,c]))
    else:
        z_1 = np.max(np.roots([a,b,c]))
    
    r_1 = r_0 + z_1 * np.tan(theta_0)
    theta_z = -np.arcsin(r_1/R)
    theta_n = -theta_z + theta_0
    
    if D < 2*r_1:
        print('The ray is out of bounds for your optical component')

    return z_1,r_1,theta_n,theta_z



def propagate2flat(r_0,theta_0,d,D):
    """
    PROPAGATE2FLAT determines the location (z_1, r_1) of a ray with initial
    height (r_0) and angle (theta_0) that hits a flat interface of diameter
    D at a distance d. Also, the angle to the surface's normal (theta_n) and
    the angle of the normal to the optical axis are calculated.
    """
    
    r_1 = r_0 + np.tan(theta_0)*d
    z_1 = d
    theta_n = theta_0
    theta_z = 0
    
    if D < 2*r_1:
        print('The ray is out of bounds for your optical component')
        
    return z_1,r_1,theta_n,theta_z 

def calc_n_K(l_nm):
    l = l_nm*10**(-3)
    pre_n = 1.1273555*l**2/(l**2 - 0.00720341707) + \
            0.124412303*l**2/(l**2 - 0.0269835916) + \
            0.827100531*l**2/(l**2 - 100.384588)
    n = np.sqrt(pre_n + 1)
    return n

def calc_n_F(l_nm):
    l = l_nm*10**(-3)
    pre_n = 1.34533359*l**2/(l**2 - 0.00997743871) + \
            0.209073176*l**2/(l**2 - 0.0470450767) + \
            0.937357162*l**2/(l**2 - 111.886764)
    n = np.sqrt(pre_n + 1)
    return n

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B)