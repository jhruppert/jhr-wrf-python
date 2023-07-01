#!/usr/bin/env python
# coding: utf-8

# Miscellaneous thermodynamic functions.
# 
# Most of the moisture calculations lean on IDL code written by
#   Dominik Brunner (brunner@atmos.umnw.ethz.ch), August 2001,
#   converted to python here.
# 
# James Ruppert  
# jruppert@ou.edu  
# May 2022


import numpy as np


############################################################################
## BASIC VARIABLES #########################################################
############################################################################

## Potential temp ######################################################

# Calculate potential temperature
#   tmpk - temp [K]
#   pres - pressure [Pa]
def theta_dry(T, pres):
    
    if np.max(pres) < 1e4:
        pres*=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    p0=1.e5 # Pa
    rd=287.04 # J/K/kg
    cp=1004. # J/K/kg
    rocp = rd/cp
    return T * ( p0 / pres ) ** rocp


## Virtual potential temp ######################################################

# Calculate virtual potential temperature
#   tmpk - temp [K]
#   qv   - vapor mixing ratio [kg/kg]
#   pres - pressure [Pa]
def theta_virtual(T, qv, pres):
    
    if np.max(pres) < 1e4:
        pres*=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    p0=1.e5 # Pa
    rd=287.04 # J/K/kg
    cp=1004. # J/K/kg
    rocp = rd/cp
    virt_corr = (1. + 0.61*qv)
    return T * virt_corr * ( p0 / pres ) ** rocp


## Density moist ######################################################

# Calculate density for an array in pressure coordinates
#   tmpk - temp [K]
#   qv   - water vapor mixing ratio [kg/kg]
#   pres - pressure [Pa]
def density_moist(T, qv, pres):
    
    if np.max(pres) < 1e4:
        pres*=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    rd=287.04
    # rv=461.5
    # eps_r=rv/rd
    # return pres / ( rd * T * (1. + qv*eps_r)/(1.+qv) )
    return pres / ( rd * T * (1. + 0.61*qv) )


## Density dry ######################################################

# Calculate density for an array in pressure coordinates
#   tmpk - temp [K]
#   pres - pressure [Pa]
def density_dry(T, pres):
    
    if np.max(pres) < 1e4:
        pres*=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    rd=287.04
    return pres / ( rd * T )


############################################################################
## MOISTURE VARIABLES ######################################################
############################################################################


## Equivalent potential temperature ##############################################

# ; BASED ON (34) OF BRYAN & FRITSCH 2002, OR (2.31) OF MARKOWSKI AND RICHARDSON
# ; (2002), which is the "wet equivalent potential temperature" (BF02) or simply
# ; "equiv pot temp" (MR02).
# ;
# ; INPUTS:
#     T - temp [K]
#     rv   - water vapor mixing ratio [kg/kg]
#     pres - pressure [Pa]
# ; XX  RTOT: TOTAL WATER (VAPOR+HYDROMETEOR) MIXING RATIO (KG/KG)
# ; 
# ; RETURNS:
# ; 
# ;   EQUIVALENT OTENTIAL TEMPERATURE (K)
# ; 
# ; James Ruppert, jruppert@ou.edu
# ; 8/4/14
# ; Converted to python, June 2022
# 
def theta_equiv(T, rv, rtot, pres):
    
    if np.max(pres) < 1e4:
        pres*=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
  # ;CONSTANTS
    R=287.    # J/K/kg
    lv0=2.5e6 # J/kg
    cp=1004.  # J/K/kg
    cpl=4186. # J/k/kg
    cpv=1885. # J/K/kg
    eps=18.0160/28.9660 # Mw / Md (source: Brunner scripts)

  # ;LATENT HEAT OF VAPORIZATION
    lv = lv0 - (cpl-cpv)*(T-273.15)

  # ;DRY AIR PRESSURE
    e = pres / ((eps/rv) + 1.)
    p_d = pres-e

  # ;CALCULATE THETA-E
    c_term = cp + cpl*rtot
    th_e = T * (1e5/p_d)**(R/c_term) * np.exp( lv*rv / (c_term*T) )

    return th_e


## Mixing ratio from vapor pressure ######################################################

# ; PURPOSE:
# ;       Convert vapor pressure (e; Pa) to mixing ratio (kg/kg).
# ; INPUTS:
# ;       e: Float or FltArr(n) vapor pressure (Pa)
# ;       p: Float or FltArr(n) ambient pressure (Pa)
# ; OUTPUTS:
# ;       rv: mixing ratio (kg/kg)
# ;                                      Mw*e              e
# ;  W (mixing ratio) = m_h2o/m_dry = -------- = Mw/Md * ---
# ;                                    Md*(p-e)           p-e
#   James Ruppert
#   8 Jan 2023
#   jruppert@ou.edu

def mixr_from_e(e,p):

    Mw=18.0160 # molecular weight of water
    Md=28.9660 # molecular weight of dry air

    rv = Mw/Md * e / (p - e) # kg/kg

    return(rv)


## Relative humidity (including for ice) ######################################################

# ; PURPOSE:
# ;       Convert mixing ratio (kg H2O per kg of dry air) at given
# ;       temperature and pressure into relative humidity (%)
# ; INPUTS:
# ;       MIXR: Float or FltArr(n) H2O mixing ratios in kg H2O per kg dry air
# ;       p   : Float or FltArr(n) ambient pressure in Pa
# ;       T   : Float or FltArr(n) ambient Temperature in C or K
# ;     ice   : switch (set = 1 to turn on) to calculate saturation over ice where T < 273.16 K
# ; OUTPUTS:
# ;       returns the relative humidity over liquid water or over ice
# ;       (if keyword /ice is set)
# ; MODIFICATION HISTORY:
# ;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), August 2001

# James Ruppert (jruppert@ou.edu), converted to python and placed here, June 2022
#   Converted all input/output to SI units, June 2022
#   Added switch and if-statement for ice, June 2022

# ;  Derivation:
# ;                                      Mw*e              e
# ;  W (mixing ratio) = m_h2o/m_dry = -------- = Mw/Md * ---
# ;                                    Md*(p-e)           p-e
# ;
# ;  RH (rel. hum.)    = e/esat(T)*100.

def relh(MIXR,p,T,ice):
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    es=esat(T)
    # if ice == 1:
    #     es[(T < 273.16)]=eice(T[(T < 273.16)])
    
    Mw=18.0160 # molecular weight of water
    Md=28.9660 # molecular weight of dry air
    fact=MIXR*Md/Mw
    return p/es*fact/(1+fact)*100.


## Saturation vapor pressure ######################################################

# ; PURPOSE:
# ;       compute saturation vapor pressure given temperature in K or C
# ; INPUTS:
# ;       T       SCALAR OR VECTOR OF TEMPERATURES IN CELSIUS OR K
# ; OUTPUTS:
# ;       returns the saturation vapor pressure in Pa
# ; MODIFICATION HISTORY:
# ;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), Feb 2000
# ;       A good reference is Gibbins, C.J., Ann. Geophys., 8, 859-886, 1990

# James Ruppert (jruppert@ou.edu), converted to python and placed here, June 2022
#   Converted all input/output to SI units, June 2022

# ; Formula with T = temperature in K
# ;    esat = exp( -6763.6/(T+T0) - 4.9283*alog((T+T0)) + 54.2190 )
# ; Formula close to that of Magnus, 1844 with temperature TC in Celsius
# ;    ESAT = 6.1078 * EXP( 17.2693882 * TC / (TC + 237.3) ) ; TC in Celsius
# ; or Emanuel's formula (also approximation in form of Magnus' formula,
# ; 1844), which was taken from Bolton, Mon. Wea. Rev. 108, 1046-1053, 1980.
# ; This formula is very close to Goff and Gratch with differences of
# ; less than 0.25% between -50 and 0 deg C (and only 0.4% at -60degC)    
# ;    esat=6.112*EXP(17.67*TC/(243.5+TC))
# 
# ; WMO reference formula is that of Goff and Gratch (1946), slightly
# ; modified by Goff in 1965:

def esat(T):
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    e1=101325.0
    TK=273.16
    # esat=e1*10**(10.79586*(1-TK/T)-5.02808*np.log10(T/TK)+
    #             1.50474*(1e-4)*(1-10**(-8.29692*(T/TK-1)))+
    #             0.42873*(1e-3)*(10**(4.76955*(1-TK/T))-1)-2.2195983)


    # JHR (11/2022): Revised based on http://cires1.colorado.edu/~voemel/vp.html
    # esat=e1*10**(10.79574*(1-TK/T)
    #              -5.02800*np.log10(T/TK)
    #              +1.50475*(1e-4)*(1-(10**(-8.2969*(T/TK-1))))
    #              +0.42873*(1e-3)*(10**(4.76955*(1-TK/T))-1)
    #              +0.78614)

    # Bolton 1980 (much faster)
    # need T in Celsius
    # esat=611.2*np.exp(17.67*(T-TK)/(243.5+(T-TK)))

    # WMO, 2008
    # http://cires1.colorado.edu/~voemel/vp.html
    # https://library.wmo.int/doc_num.php?explnum_id=7450
    esat=611.2*np.exp(17.62*(T-TK)/(243.12+(T-TK)))

    return esat


## Saturation vapor pressure over ice ######################################################

# ; PURPOSE:
# ;       compute saturation vapor pressure over ice given temperature
# ;       in K or C. The accuracy is stated by Marti & Mauersberger
# ;       as 2% in the temperature range of 170K to 250K.
# ; INPUTS:
# ;       T       SCALAR OR VECTOR OF TEMPERATURES IN CELSIUS OR K
# ; OUTPUTS:
# ;       returns the saturation vapor pressure in Pa
# ; MODIFICATION HISTORY:
# ;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), March 2000
# ;       reference: Marti and Mauersberger, GRL 20, 363-366, 1993.

# James Ruppert (jruppert@ou.edu), converted to python and placed here, June 2022
#   Converted all input/output to SI units, June 2022

def eice(T):
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    # ; Define constants
    # A=-2663.5
    # B=12.537
    # logp=A/T+B
    # # return (10**logp)/100. # conversion to hPa
    # return 10**logp
    
    # Guide to Meteorological Instruments and Methods of Observation (CIMO Guide)
    # -- WMO, 2008
    # http://cires1.colorado.edu/~voemel/vp.html
    # https://library.wmo.int/doc_num.php?explnum_id=7450
    TK=273.16
    eice=611.2*np.exp(22.46*(T-TK)/(272.62+(T-TK)))
    return eice



