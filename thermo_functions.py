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
from metpy.calc import cape_cin, dewpoint_from_specific_humidity, parcel_profile#, most_unstable_parcel
from metpy.units import units


############################################################################
## BASIC VARIABLES #########################################################
############################################################################

## Potential temp ######################################################

# Calculate potential temperature
#   tmpk - temp [K]
#   pres - pressure [Pa]
def theta_dry(T, pres):

    p_fact=1
    if np.max(pres) < 1e4:
        p_fact=1e2 # Convert to Pa

    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0

    p0=1e5 # Pa
    rd=287.04 # J/K/kg
    cp=1004. # J/K/kg
    rocp = rd/cp
    return T * ( p0 / (pres*p_fact) ) ** rocp


## Virtual potential temp ######################################################

# Calculate virtual potential temperature
#   tmpk - temp [K]
#   qv   - vapor mixing ratio [kg/kg]
#   pres - pressure [Pa]
def theta_virtual(T, qv, pres):

    p_fact=1
    if np.max(pres) < 1e4:
        p_fact=1e2 # Convert to Pa

    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0

    p0=1e5 # Pa
    rd=287.04 # J/K/kg
    cp=1004. # J/K/kg
    rocp = rd/cp
    # virt_corr = (1. + 0.61*qv)
    return T * (1. + 0.61*qv) * ( p0 / (pres*p_fact) ) ** rocp


## Density moist ######################################################

# Calculate density for an array in pressure coordinates
#   T    - temp [C or K]
#   qv   - water vapor mixing ratio [kg/kg]
#   pres - pressure [Pa]
def density_moist(T, qv, pres):

    p_fact=1
    if np.max(pres) < 1e4:
        p_fact=1e2 # Convert to Pa

    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0

    rd=287.04
    # rv=461.5
    # eps_r=rv/rd
    # return pres / ( rd * T * (1. + qv*eps_r)/(1.+qv) )
    return pres*p_fact / ( rd * T * (1. + 0.61*qv) )


## Density dry ######################################################

# Calculate density for an array in pressure coordinates
#   T    - temp [C or K]
#   pres - pressure [Pa]
def density_dry(T, pres):

    p_fact=1
    if np.max(pres) < 1e4:
        p_fact=1e2 # Convert to Pa

    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.

    rd=287.04
    return pres*p_fact / ( rd * (T+T0) )


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
    
    p_fact=1
    if np.max(pres) < 1e4:
        p_fact=1e2 # Convert to Pa
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    # T+=T0
    
  # ;CONSTANTS
    rd=287.    # J/K/kg
    lv0=2.5e6 # J/kg
    cp=1004.  # J/K/kg
    cpl=4186. # J/k/kg
    cpv=1885. # J/K/kg
    eps=18.0160/28.9660 # Mw / Md (source: Brunner scripts)

  # ;LATENT HEAT OF VAPORIZATION
    # lv = lv0 - (cpl-cpv)*(T-273.15)

  # ;DRY AIR PRESSURE
    # e = pres / ((eps/rv) + 1.)
    # p_d = pres - e
    # p_d = (pres*p_fact) - ((pres*p_fact) / ((eps/rv) + 1.))

  # ;CALCULATE THETA-E
    p0=1e5 # Pa
    # c_term = cp + cpl*rtot
    th_e = (T+T0) * (p0/ ((pres*p_fact) - ((pres*p_fact) / ((eps/rv) + 1.))) )**(rd/(cp + cpl*rtot)) \
        * np.exp( (lv0 - (cpl-cpv)*((T+T0)-273.15))*rv / ((cp + cpl*rtot)*(T+T0)) )

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

def calc_relh(MIXR,p,T,ice=True):
    
    if np.min(T) < 105.: # degC or K?
        T0=273.16
    else:
        T0=0.
    T+=T0
    
    es=esat(T)
    if ice:
        es[(T < 273.16)]=eice(T[(T < 273.16)])
    
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

    # e1=101325.0
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
    # Put T into K
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

    # Trying Buck's formula (1996)
    # https://www.eas.ualberta.ca/jdwilson/EAS372_13/Vomel_CIRES_satvpformulae.html
    # http://cires1.colorado.edu/~voemel/vp.html
    # eice=611.15*np.exp((23.036 - (T-TK)/333.7)*(T-TK)/(279.82+(T-TK)))
    return eice


## Saturation mixing ratio ######################################################

# Just combining a couple chunks of code from other functions in this routine.
# 
# All input/output assumed SI units, June 2022
# 
# James Ruppert (jruppert@ou.edu)

def rv_saturation(T, pres):

    if np.min(T) < 105.: # degC or K?
        T0=0
    else:
        T0=-273.16

    esat=611.2*np.exp(17.62*(T+T0)/(243.12+(T+T0))) # Pa

    Mw=18.0160 # molecular weight of water
    Md=28.9660 # molecular weight of dry air

    rv_sat = Mw/Md * esat / (pres - esat) # kg/kg

    return rv_sat

## Specific humidity from mixing ratio ######################################################

# ; PURPOSE:
# ;       Convert mixing ratio (kg H2O per kg of dry air) to
# ;       specific humidity (same units).
# ; INPUTS:
# ;       MIXR: Float or FltArr(n) H2O mixing ratios in kg H2O per kg dry air
# ; OUTPUTS:
# ;       returns the specific humidity
# 
# James Ruppert (jruppert@ou.edu), September 2024, from the mid-Atlantic

def mixr2sh(mixr):
    q = mixr / (1 + mixr)
    return q


## CAPE / CIN ######################################################

# Use MetPy package to calculate CAPE and CIN from the most unstable parcel.
# https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.cape_cin.html
# 
# Uses MetPy packages:
#   dewpoint_from_specific_humidity
#   parcel_profile
#   cape_cin
# 
# Inputs:
#   tmpk - temp in [ºC or K]
#   qv   - water vapor mixing ratio [kg/kg]
#   pres - pressure [hPa or Pa]
# 
# Options:
#   type = 'mu' or 'sfc' for 1) most unstable or 2) surface parcel
#       default: 'sfc'
# 
# Dimensions of input:
#   Either (nz) or (nt, nz) where nt refers to time.
# 
# Returns:
#   CAPE - J/kg
#   CIN  - J/kg
# 
# James Ruppert
# jruppert@ou.edu
# September 2024 from the Meteor
#  
def get_cape_cin(T_in, qv_in, pres_in, type='sfc'):

    # Convert hPa if necessary
    p_fact=1
    if np.nanmax(pres_in) > 1e5:
        p_fact=1e-2 # Convert to hPa

    # Convert T to ºC if necessary
    if np.nanmin(T_in) < 105.: # degC or K?
        T0=0
    else:
        T0=273.16

    # Check for time dimension
    nt = 1
    if T_in.ndim > 1:
        nt = T_in.shape[0]

    cape = np.zeros(nt)
    cin  = np.zeros(nt)

    # Get specific humidity from mixing ratio
    sh_in = mixr2sh(qv_in)
    # Get theta_e for finding most unstable parcel
    theta_e_in = theta_equiv(T_in, qv_in, qv_in, pres_in)

    for it in range(nt):

        if nt > 1:
            pres = pres_in[it]
            T = T_in[it]
            sh = sh_in[it]
            theta_e = theta_e_in[it]
        else:
            pres = pres_in
            T = T_in
            sh = sh_in
            theta_e = theta_e_in

        # Check for, remove NaNs
        ind = np.isfinite(pres)
        # Skip entire sounding if all NaN
        if np.count_nonzero(ind) == 0:
            cape[it] = np.nan
            cin[it] = np.nan
            continue
        p = pres[ind] * p_fact
        tmpc = T[ind] - T0
        sh = sh[ind]
        theta_e = theta_e[ind]

        # Check if masked array
        if isinstance(tmpc, np.ma.MaskedArray):
            p    = np.ma.getdata(p)
            tmpc = np.ma.getdata(tmpc)
            sh   = np.ma.getdata(sh)
            theta_e = np.ma.getdata(theta_e)

        # Apply units for MetPy functions
        p *= units.hPa
        tmpc *= units.degC
        sh *= units.dimensionless

        # Run MetPy routines
        dwpc = dewpoint_from_specific_humidity(p, tmpc, sh)
        # mu_parcel = most_unstable_parcel(p, tmpc, dwpc, depth=50*units.hPa)

        if type == 'mu':
            # Use maximum* theta-e as most unstable parcel
            #  * checked within 50 hPa of surface pressure
            p_check = np.array(p)
            iz_check = np.where(p_check >= np.max(p_check)-50)[0]
            iz_parcel = np.where(theta_e[iz_check] == np.max(theta_e[iz_check]))[0][0]
        else:
            iz_parcel=0

        # In case of no unstable parcel found
        try:
            mu_prof = parcel_profile(p, tmpc[iz_parcel], dwpc[iz_parcel]).to('degC')
            icape, icin = cape_cin(p, tmpc, dwpc, mu_prof, which_lfc='bottom', which_el='top')
            cape[it] = np.array(icape)
            cin[it] = np.array(icin)
        except:
            cape[it] = np.nan
            cin[it] = np.nan

    return cape, cin