#!/usr/bin/env python
# coding: utf-8

# ### Miscellaneous functions (mostly thermodynamic)
# 
# James Ruppert  
# jruppert@ou.edu  
# 5/27/22

import numpy as np

## Potential temp ######################################################

# Calculate potential temperature
#   tmpk   - temp in [K]
#   pres1d - pressure as 1D array in [Pa]
def theta_dry(tmpk, pres1d):
    p0=1.e5 # Pa
    rd=287.04 # J/K/kg
    cp=1004. # J/K/kg
    rocp = rd/cp
    return tmpk * ( p0 / pres1d[np.newaxis,:,np.newaxis,np.newaxis] ) ** rocp

## Density moist ######################################################

# Calculate density for an array in pressure coordinates
#   Assumes pres is 1D and other vars are 4d, with vertical in the 2nd dimension
#   tmpk   - temp in [K]
#   qv     - water vapor mixing ratio [kg/kg]
#   pres1d - pressure as 1D array in [Pa]
def density_moist(tmpk, qv, pres1d):
    rd=287.04
    rv=461.5
    eps_r=rv/rd
    # virt_corr = (1. + qv*eps_r)/(1.+qv)
    virt_corr = (1. + 0.61*qv)
    return pres1d[np.newaxis,:,np.newaxis,np.newaxis] / ( rd * tmpk * virt_corr )

## Density dry ######################################################

# Calculate density for an array in pressure coordinates
#   Assumes pres is 1D and other vars are 4d, with vertical in the 2nd dimension
#   tmpk   - temp in [K]
#   pres1d - pressure as 1D array in [Pa]
def density_dry(tmpk, pres1d):
    rd=287.04
    return pres1d[np.newaxis,:,np.newaxis,np.newaxis] / ( rd * tmpk )
