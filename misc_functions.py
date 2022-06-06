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


## Relative humidity ######################################################

# Calculate RH wrt ice where T < 0C
#   Assumes pres is 1D and other vars are 4d, with vertical in the 2nd dimension
#   tmpk   - temp in [K]
#   qv     - water vapor mixing ratio [kg/kg]
#   pres1d - pressure as 1D array in [Pa]
def density_dry(tmpk, qv, pres1d):
    rd=287.04
    return pres1d[np.newaxis,:,np.newaxis,np.newaxis] / ( rd * tmpk )


## Saturation vapor pressure ######################################################

# ; PURPOSE:
# ;       compute saturation vapor pressure given temperature in K or C
# ; INPUTS:
# ;       T       SCALAR OR VECTOR OF TEMPERATURES IN CELSIUS OR K
# ; OUTPUTS:
# ;       returns the saturation vapor pressure in hPa
# ; MODIFICATION HISTORY:
# ;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), Feb 2000
# ;       A good reference is Gibbins, C.J., Ann. Geophys., 8, 859-886, 1990
# ;  James Ruppert (jruppert@ou.edu), converted to python and placed here, June 2022

# ; Formula with T = temperature in K
# ;    esat = exp( -6763.6/(T+T0) - 4.9283*alog((T+T0)) + 54.2190 )
# 
# ; Formula close to that of Magnus, 1844 with temperature TC in Celsius
# ;    ESAT = 6.1078 * EXP( 17.2693882 * TC / (TC + 237.3) ) ; TC in Celsius
# 
# ; or Emanuel's formula (also approximation in form of Magnus' formula,
# ; 1844), which was taken from Bolton, Mon. Wea. Rev. 108, 1046-1053, 1980.
# ; This formula is very close to Goff and Gratch with differences of
# ; less than 0.25% between -50 and 0 deg C (and only 0.4% at -60degC)    
# ;    esat=6.112*EXP(17.67*TC/(243.5+TC))
# 
# ; WMO reference formula is that of Goff and Gratch (1946), slightly
# ; modified by Goff in 1965:

def esat(T):
  if T[0] < 105.:
     THEN T0=273.16 ELSE T0=0.
    e1=1013.250
    TK=273.16
    esat=e1*10^(10.79586*(1-TK/(T+T0))-5.02808*alog10((T+T0)/TK)+$
                1.50474*1e-4*(1-10^(-8.29692*((T+T0)/TK-1)))+$
                0.42873*1e-3*(10^(4.76955*(1-TK/(T+T0)))-1)-2.2195983)

    return esat


# Calculate esat RH wrt ice where T < 0C
#   Assumes pres is 1D and other vars are 4d, with vertical in the 2nd dimension
#   tmpk   - temp in [K]
#   qv     - water vapor mixing ratio [kg/kg]
#   pres1d - pressure as 1D array in [Pa]
def density_dry(tmpk, qv, pres1d):
    rd=287.04
    return pres1d[np.newaxis,:,np.newaxis,np.newaxis] / ( rd * tmpk )




