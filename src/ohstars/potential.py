# ----------------------------------------------------------------------------
#
# TITLE - potential.py
# AUTHOR - James Lane
# PROJECT - OHStars
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Tools for dealing with potentials.

Computation of McMillan 2017 potential from: 
https://github.com/jmackereth/orbit-estimation/
'''

__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb

## galpy
from galpy import potential
from galpy.util import bovy_conversion as gpconv

# ----------------------------------------------------------------------------

# Module globals
ro=8.21
vo=233.1
sigo = gpconv.surfdens_in_msolpc2(vo=vo, ro=ro)
rhoo = gpconv.dens_in_msolpc3(vo=vo, ro=ro)

#gas disk parameters (fixed in McMillan (2017)...)
Rd_HI = 7/ro
Rm_HI = 4/ro
zd_HI = 0.085/ro
Sigma0_HI = 53.1/sigo
Rd_H2 = 1.5/ro
Rm_H2 = 12/ro
zd_H2 = 0.045/ro
Sigma0_H2 = 2180/sigo

#parameters of best-fitting model in McMillan (2017)
#stellar disks
Sigma0_thin = 896/sigo
Rd_thin = 2.5/ro
zd_thin = 0.3/ro
Sigma0_thick = 183/sigo
Rd_thick = 3.02/ro
zd_thick = 0.9/ro
#bulge
rho0_bulge = 98.4/rhoo
r0_bulge = 0.075/ro
rcut = 2.1/ro
#DM halo
rho0_halo = 0.00854/rhoo
rh = 19.6/ro

class mySCFPotential(potential.SCFPotential):
    '''mySCFPotential:
    '''
    def _R2deriv(self,R,z,phi=0.,t=0.):
        dR= 1e-8
        return (self._Rforce(R,z) - self._Rforce(R+dR,z))/dR
    #def
    
    def _z2deriv(self,R,z,phi=0.,t=0.):
        dz = 1e-8
        return (self._zforce(R,z) - self._zforce(R,z+dz))/dz
    #def
    
    def _Rzderiv(self,R,z,phi=0.,t=0.):
        dR = 1e-8
        return (self._zforce(R,z) - self._zforce(R+dR,z))/dR
    #def
#cls

class myDiskSCFPotential(potential.DiskSCFPotential):
    '''myDiskSCFPotential:
    '''
    def _R2deriv(self,R,z,phi=0.,t=0.):
        dR= 1e-8
        return (self._Rforce(R,z) - self._Rforce(R+dR,z))/dR
    #def
    
    def _z2deriv(self,R,z,phi=0.,t=0.):
        dz = 1e-8
        return (self._zforce(R,z) - self._zforce(R,z+dz))/dz
    #def
    
    def _Rzderiv(self,R,z,phi=0.,t=0.):
        dR = 1e-8
        return (self._zforce(R,z) - self._zforce(R+dR,z))/dR
    #def
#cls

def gas_dens(R,z, separate=False):
    '''gas_dens:
    
    Calculate gas disk density as a f(R,z).
    
    May not handle R/z == 0 correctly... need to avoid divide by 0.
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        separate (bool) - Return HI and H2 density as a pair [False]
        
    Returns:
        Density (float or 2-array)
    '''
    if hasattr(R, "__iter__"):
        HI_dens = np.empty(len(R))
        H2_dens = np.empty(len(R))
        Requal0 = R == 0.
        HI_dens[Requal0] = 0.
        H2_dens[Requal0] = 0.
        HI_dens[~Requal0] = Sigma0_HI/(4*zd_HI)*np.exp(-Rm_HI/R[~Requal0]-R[~Requal0]/Rd_HI)*(sech(z[~Requal0]/(2*zd_HI)))**2
        H2_dens[~Requal0] = Sigma0_H2/(4*zd_H2)*np.exp(-Rm_H2/R[~Requal0]-R[~Requal0]/Rd_H2)*(sech(z[~Requal0]/(2*zd_H2)))**2
        if separate:
            return HI_dens, H2_dens
        ##fi
        return HI_dens+H2_dens   
    else:
        if R == 0.:
            HI_dens = Sigma0_HI/(4*zd_HI)*np.exp(-np.inf-R/Rd_HI)*(sech(z/(2*zd_HI)))**2
            H2_dens = Sigma0_H2/(4*zd_H2)*np.exp(-np.inf-R/Rd_H2)*(sech(z/(2*zd_H2)))**2
            HI_dens = 0.
            H2_dens = 0.
        else:
            HI_dens = Sigma0_HI/(4*zd_HI)*np.exp(-Rm_HI/R-R/Rd_HI)*(sech(z/(2*zd_HI)))**2
            H2_dens = Sigma0_H2/(4*zd_H2)*np.exp(-Rm_H2/R-R/Rd_H2)*(sech(z/(2*zd_H2)))**2
        ##ie
        if separate:
            return HI_dens, H2_dens
        ##fi
        return HI_dens+H2_dens
    ##ie
#def

def stellar_dens(R,z):
    '''stellar_dens:
    
    Calculate the thick & thin stellar density
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)
    '''
    thin_dens = Sigma0_thin/(2*zd_thin)*np.exp(-np.fabs(z)/zd_thin-R/Rd_thin)
    thick_dens = Sigma0_thick/(2*zd_thick)*np.exp(-np.fabs(z)/zd_thick-R/Rd_thick)
    return thin_dens+thick_dens
#def

def bulge_dens(R,z):
    '''bulge_dens:
    
    Calculate bulge density
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)
    '''
    rdash = np.sqrt(R**2+(z/0.5)**2)
    dens = rho0_bulge/(1+rdash/r0_bulge)**1.8*np.exp(-(rdash/rcut)**2)
    return dens
#def

def NFW_dens(R,z):
    '''NFW_dens:
    
    Calculate the NFW as f(R,z), not used to generate the potential
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)
    
    '''
    r = np.sqrt(R**2+z**2)
    x = r/rh
    dens = rho0_halo/(x*(1+x)**2)
    return dens
#def

def bulge_gas_dens(R,z):
    '''bulge_gas_dens:
    
    Combined bulge, gas and stellar density models
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)
    '''
    return bulge_dens(R,z)+gas_dens(R,z)+stellar_dens(R,z)
#def

def gas_stellar_dens(R,z):
    '''gas_stellar_dens:
    
    Combine only the gas and stellar densities
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)    
    '''
    return gas_dens(R,z)+stellar_dens(R,z)
#def

def tot_dens(R,z):
    '''tot_dens:
    
    Calculat the total density, including NFW halo, as f(R,z)
    
    Args:
        R (float) - Galactocentric R
        z (float) - Galactocentric z
        
    Returns:
        Density (float)
    '''
    return gas_dens(R,z)+stellar_dens(R,z)+bulge_dens(R,z)+NFW_dens(R,z)
#def

def sech(x):
    '''sech:
    
    Trig helper    
    '''
    return 1./np.cosh(x)
#def

def make_McMillan2017(ro=8.21, vo=233.1):
    '''
    make_McMillan2017:
    
    Make the McMillan 2017 potential using SCF basis expansion:
    
    No arguments, but could change ro, vo if desired. ro is from the Gravity 
    Collaboration (2018). vo is from Eilers (2018) and Schonrich (2010)
    
    Args:
        None
    '''

    #dicts used in DiskSCFPotential 
    sigmadict = [{'type':'exp','h':Rd_HI,'amp':Sigma0_HI, 'Rhole':Rm_HI},
                 {'type':'exp','h':Rd_H2,'amp':Sigma0_H2, 'Rhole':Rm_H2},
                 {'type':'exp','h':Rd_thin,'amp':Sigma0_thin, 'Rhole':0.},
                 {'type':'exp','h':Rd_thick,'amp':Sigma0_thick, 'Rhole':0.}]

    hzdict = [{'type':'sech2', 'h':zd_HI},
              {'type':'sech2', 'h':zd_H2},
              {'type':'exp', 'h':0.3/ro},
              {'type':'exp', 'h':0.9/ro}]
    McMillan_bulge=\
    mySCFPotential(Acos=potential.scf_compute_coeffs_axi(bulge_dens,20,10,a=0.1)[0],
                    a=0.1,ro=ro,vo=vo)
    McMillan_disk = myDiskSCFPotential(dens=lambda R,z: gas_stellar_dens(R,z),
                                     Sigma=sigmadict, hz=hzdict,
                                     a=2.5, N=30, L=30,ro=ro,vo=vo)
    McMillan_halo = potential.NFWPotential(amp = rho0_halo*(4*np.pi*rh**3),
                                 a = rh,ro=ro,vo=vo)
    McMillan2017 = [McMillan_disk,McMillan_halo,McMillan_bulge]
    
    return McMillan2017
#def

# Make the potential so it's just an attribute.

#dicts used in DiskSCFPotential 
sigmadict = [{'type':'exp','h':Rd_HI,'amp':Sigma0_HI, 'Rhole':Rm_HI},
             {'type':'exp','h':Rd_H2,'amp':Sigma0_H2, 'Rhole':Rm_H2},
             {'type':'exp','h':Rd_thin,'amp':Sigma0_thin, 'Rhole':0.},
             {'type':'exp','h':Rd_thick,'amp':Sigma0_thick, 'Rhole':0.}]

hzdict = [{'type':'sech2', 'h':zd_HI},
          {'type':'sech2', 'h':zd_H2},
          {'type':'exp', 'h':0.3/ro},
          {'type':'exp', 'h':0.9/ro}]
McMillan_bulge=\
mySCFPotential(Acos=potential.scf_compute_coeffs_axi(bulge_dens,20,10,a=0.1)[0],
                a=0.1,ro=ro,vo=vo)
McMillan_disk = myDiskSCFPotential(dens=lambda R,z: gas_stellar_dens(R,z),
                                 Sigma=sigmadict, hz=hzdict,
                                 a=2.5, N=30, L=30,ro=ro,vo=vo)
McMillan_halo = potential.NFWPotential(amp = rho0_halo*(4*np.pi*rh**3),
                             a = rh,ro=ro,vo=vo)
McMillan2017 = [McMillan_disk,McMillan_halo,McMillan_bulge]

# ----------------------------------------------------------------------------
