# ----------------------------------------------------------------------------
#
# TITLE - analysis
# AUTHOR - James Lane
# PROJECT - OHStars
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Functions for analysis
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb
# import copy
# import glob
# import subprocess

## Plotting
from matplotlib import pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import colors
# from matplotlib import cm
# import aplpy

## Astropy
# from astropy.io import fits
# from astropy.coordinates import SkyCoord
# from astropy import table
# from astropy import units as u
# from astropy import wcs

## galpy
# from galpy import orbit
# from galpy import potential
# from galpy.util import bovy_coords as gpcoords
# from galpy.util import bovy_conversion as gpconv
# from galpy.util import bovy_plot as gpplot

# ----------------------------------------------------------------------------

def sample_kinematics( data,
                       n_samples = 1000,
                    ):
    '''
    sample_kinematics
    
    Sample the 5D kinematics from gaia data for a single 
    star. Takes into account all covariances
    
    Args:
        data (astropy row) - A queriable Astropy row object with gaia_source 
            keywords
        n_samples (int) - The number of samples to draw [1000]
    
    Returns:
        samples (mxn array) - The m samples by n variables array
        cov (5x5 array) - The covariance matrix used
    '''
    
    # Assume the structure of the gaia_source data
    ps_abbr = np.array(['ra','dec','parallax','pmra','pmdec'])
    
    # Gather the means
    means = np.zeros(5)
    for j in range(5):
        if j == 0 or j == 1:
            means[j] = data[ps_abbr[j]] * (1000*3600) # Deg -> mas
        else:
            means[j] = data[ps_abbr[j]]
        ##fi
    ###j
    
    # Generate the covaraiance matrix
    cov = np.zeros((5,5))
    for j in range(5):
        for k in range(5):
            if j == k:
                cov[j,k] = data[ps_abbr[j]+'_error']**2
            else:
                try:
                    corr = data[ps_abbr[j]+'_'+ps_abbr[k]+'_corr']
                except KeyError:
                    corr = data[ps_abbr[k]+'_'+ps_abbr[j]+'_corr']
                cov[j,k] = data[ps_abbr[j]+'_error'] * data[ps_abbr[k]+'_error'] * corr
            ##ie
        ###k
    ###j
    
    # Sample from the covariance matrix
    samples = np.random.multivariate_normal( means, cov, size=n_samples )
    sample_labels = np.array([r'$\alpha_{off}$ [mas]',
                              r'$\delta_{off}$ [mas]',
                              r'$\pi$',
                              r'$\mu_{\alpha}$ [mas/yr]',
                              r'$\mu_{\delta}$ [mas/yr]'])
                              
    # Transform RA/Dec from mas to degrees
    samples[:,0] /= (1000*3600)
    samples[:,1] /= (1000*3600)
    
    return samples, cov
#def

# ----------------------------------------------------------------------------
