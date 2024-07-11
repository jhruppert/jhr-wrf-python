#!/usr/bin/env python
# coding: utf-8

# Python script to calculate multiple diagnostics from WRF output, including
# MSE budget terms, and write them out.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/24/23

from netCDF4 import Dataset
import numpy as np
import os
from thermo_functions import *
from read_functions import *
from write_ncfile import write_ncfile
from time import time as runtimer
import sys


#### Main settings

# Check for required number of processors
# if nproc != nvars:
#     print("Check NPROC (-n option)! Should be ",nvars)
#     print("Killing batch job")
#     sys.exit()

storm = 'haiyan'
storm = 'maria'

filename_out='mse_diag.nc' # this is for ALL variables in the var_names list

msetop = 100 # top for MSE integrals

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    # tests = ['ctl','ncrf36h']
    # tests = ['crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['STRATANVIL_OFF','STRAT_OFF']
    tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    tests = ['ctl','ncrf48h']#'ncrf36h']
    tests = [tests[1],'crfon72h']
    tests = ['crfon72h']
    tests = ['ctl','ncrf48h','ncrf36h']
ntest=len(tests)

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 20
# nmem = 1

######################################################################

# Ens member strings
memb0=1 # Starting member to read
nums=np.arange(memb0,nmem+memb0,1)
nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

datdir2 = 'post/d02/'

##### Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
datdir3d = datdir #+'v2/'
nt, nz, nx1, nx2, pres = get_file_dims(datdir)
dp = (pres[0]-pres[1])*1e2 # Pa
kmsetop = np.where(pres == msetop)[0][0]

# Get Lat/Lon
wrfdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
dirlist = os.listdir(wrfdir)
subs="wrfout_d02"
wrf_files = list(filter(lambda x: subs in x, dirlist))
wrf_files = [wrfdir + s for s in wrf_files]
varfil_main = Dataset(wrf_files[0])
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()
lon1d=lon[0,:]
lat1d=lat[:,0]

deg2m = np.pi*6371*1e3/180
x1d = lon1d * deg2m
y1d = lat1d * deg2m

######################################################################

#### NetCDF variable read functions

def var_read(datdir,varname,ikread):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,0:ikread+1,:,:]
    varfil_main.close()
    return var

######################################################################

#### NetCDF variable metadata

def var_ncdf_metadata(dims3d):
    
    var_names = [
        'pres',
        'rho',
        'theta_e',
        'pw',
        'pw_sat',
        'vmf',
        'vmfu',
        'vmfd',
        'condh',
        'dse',
        'mse',
        'mse_vint',
        'vadv_dse_vint',
        'vadv_mse_vint',
        'hadv_dse_vint',
        'hadv_mse_vint',
        'dse_diverg_vint',
        'mse_diverg_vint',
        'dse_fluxdiverg_vint',
        'mse_fluxdiverg_vint',
    ]
    descriptions = [
        'pressure levels',
        'density of moist air',
        'equivalent potential temperature',
        'water vapor integrated over full column from postproc data',
        'saturation water vapor integrated over full column from postproc data',
        'vertical mass flux vertically integrated (up to 100 hPa)',
        'upward-masked mass flux vertically integrated (up to 100 hPa)',
        'downward-masked mass flux vertically integrated (up to 100 hPa)',
        'condensation heating from H_DIABATIC vertically int (up to 100 hPa), converted to rainfall units',
        'dry static energy, calculated as cpT + gz',
        'moist static energy, calculated as cpT + gz + L_v*q',
        'vertically int moist static energy, calculated as 1/g*integral(mse)dp up to 100 hPa',
        'vertically int VADV of DSE (up to 100 hPa)',
        'vertically int VADV of MSE (up to 100 hPa)',
        'vertically int HADV of DSE (up to 100 hPa)',
        'vertically int HADV of MSE (up to 100 hPa)',
        'vertically int DSE*(del.V) (up to 100 hPa)',
        'vertically int MSE*(del.V) (up to 100 hPa)',
        'vertically int del.(DSE*V) (up to 100 hPa)',
        'vertically int del.(MSE*V) (up to 100 hPa)',
    ]
    units = [
        'hPa',
        'kg/m^3',
        'K',
        'mm',
        'mm',
        'kg/m/s',
        'kg/m/s',
        'kg/m/s',
        'mm/day',
        'J/kg',
        'J/kg',
        'J/m^2',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
    ]
    nt,nz,nx1,nx2 = dims3d
    dims2d = (nt,nx1,nx2)
    dims2d_names = ('nt','nx1','nx2')
    dims3d_names = ('nt','nz','nx1','nx2')
    dims_set = [
        [('nz',),(nz,)],
        [dims3d_names,dims3d],
        [dims3d_names,dims3d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims3d_names,dims3d],
        [dims3d_names,dims3d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
        [dims2d_names,dims2d],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set #, dim_names

##### Various diagnostics ################################

def get_3D_diagnostics(datdir, pres, dp, nz, z_b, g):

    # Constants
    cp=1004.  # J/K/kg
    cpl=4186. # J/k/kg
    cpv=1885. # J/K/kg
    lv0=2.5e6 # J/kg

    # Water vapor
    varname='QVAPOR'
    qv = var_read(datdir,varname,nz) # kg/kg
    # nt,nz,nx1,nx2 = qv.shape
    # Temperature
    varname='T'
    tmpk = var_read(datdir,varname,nz) # K

    # Height
    varname='Z'
    z = var_read(datdir,varname,nz) # m2/s2
    z += z_b
    # Dry static energy (DSE)
    dse = cp*tmpk + g*z # J/kg

    # Density
    rho = density_moist(tmpk,qv,pres[np.newaxis,:,np.newaxis,np.newaxis,]*1e2) # kg/m3
    # Theta-e
    theta_e = theta_equiv(tmpk,qv,qv,pres[np.newaxis,:,np.newaxis,np.newaxis,]*1e2) # kg/m3

    # Calculate PW and PW_SAT for saturation fraction
    qv_sat = rv_saturation(tmpk, (pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/kg
    # Integrate over full column
    pw = vert_int(qv, dp, g) # mm
    pw_sat = vert_int(qv_sat, dp, g) # mm

    # Latent heat of vaporization
    lv = lv0 - (cpl-cpv)*(tmpk-273.15)

    # Moist static energy (MSE)
    g=9.81 # m/s2
    mse = dse + lv*qv # J/kg
    mse_vint = vert_int(mse[:,0:kmsetop+1,:,:], dp, g) # J/m^2

    return rho, theta_e, dse, mse, pw, pw_sat, mse_vint

def get_condh(datdir, dp, nz, g):
    # Microphysics (MP) latent heating
    lh = var_read(datdir,'H_DIABATIC',nz) # K/s
    condh = vert_int(lh[:,0:kmsetop,:,:], dp, g) # kg*K/(m2*s) [needs *cp to get W/m2]
    return condh

# Vertical mass flux
def get_vmf_updown(w, dp, g):
    # Mask for Up/Dn
    wu = np.ma.masked_where((w < 0), w, copy=True)
    wd = np.ma.masked_where((w > 0), w, copy=True)
    vmf = vert_int(w[:,0:kmsetop,:,:], dp, g) # kg/m/s
    vmfu = vert_int(wu[:,0:kmsetop,:,:], dp, g) # kg/m/s
    vmfd = vert_int(wd[:,0:kmsetop,:,:], dp, g) # kg/m/s
    return vmf, vmfu, vmfd

##### MSE & DSE functions ################################

# All calculations assume:
#   - Lagrangian convention, i.e., partial()/partial(t) + HADV + VADV = ...
#   - dp is a constant
#   - < > is vertical integral
# 
# Accd. to Raymond et al. (2009), Inoue and Back (2015)
# for vertically integrated terms,
#   Vert-advection = < omeg . d(var)/dp >
#   H-advection = < V . Del(var) > = < u.d(var)/dx > + < v.d(var)/dy >
#   Their sum = Total flux diverg = Del . <var*V>
#                                 = d<var*u>/dx - d<var*v>/dy

def vert_int(invar, dp, g):
    # Vertically integrate: 1/g * SUM(var*dp)
    # Negative is absorbed by dp>0
    var_sum = np.sum(invar, axis=1)*dp/g
    return var_sum

def vadv_vint(w, rho, invar, dp, g):
    omeg = w * (-1)*g*rho # omega from hydrostatic approximation
    vadv = omeg * np.gradient(invar, (-dp), axis=1)*(-1)
    vadv_vint = vert_int(vadv, dp, g)
    return vadv_vint

def hadv_vint(u, v, x1d, y1d, invar, dp, g):
    adv_x = u * np.gradient(invar, x1d, axis=3)
    adv_y = v * np.gradient(invar, y1d, axis=2)
    hadv_vint = vert_int((adv_x + adv_y), dp, g)
    return hadv_vint

def diverg_vint(u, v, x1d, y1d, invar, dp, g):
    # Divergence of phi = < phi del . V >
    dudx = np.gradient(u, x1d, axis=3) # /s
    dvdy = np.gradient(v, y1d, axis=2) # /s
    div = dudx + dvdy
    diverg = vert_int((invar * div), dp, g)
    return diverg

def tot_fdiverg_vint(u, v, x1d, y1d, invar, dp, g):
    # Total flux divergence of phi:
    #   flux_diverg = < del . phi*V >
    #       = < d/dx (phi*u) + d/dy (phi*v) >
    diverg_x = np.gradient((u * invar), x1d, axis=3)
    diverg_y = np.gradient((v * invar), y1d, axis=2)
    flux_diverg = diverg_x + diverg_y
    flux_diverg_vint = vert_int(flux_diverg, dp, g)
    return flux_diverg_vint

##### Main function ######################################

g=9.81 # m/s2

# Read base-state height once per storm
varname='ZB'
z_b = var_read(datdir,varname,nz) # m

def main_routine(datdir):

    # Get dimensions
    nt, nz, nx1, nx2, pres = get_file_dims(datdir)

    # Required variables

    rho, theta_e, dse, mse, pw, pw_sat, mse_vint = get_3D_diagnostics(datdir, pres, dp, nz, z_b, g)

    condh = get_condh(datdir, dp, nz, g) # kg*K/(m2*s) [needs *cp to get W/m2]

    # Winds
    w = var_read(datdir,'W',nz) # m/s
    u = var_read(datdir,'U',nz)
    v = var_read(datdir,'V',nz)

    vmf, vmfu, vmfd = get_vmf_updown(w, dp, g)

    # Calculate advection/divergence terms and vertically integrate

    # Vertical advection
    vadv_dse_vint = vadv_vint(w[:,0:kmsetop+1,:,:], rho[:,0:kmsetop+1,:,:],
                            dse[:,0:kmsetop+1,:,:], dp, g)
    vadv_mse_vint = vadv_vint(w[:,0:kmsetop+1,:,:], rho[:,0:kmsetop+1,:,:],
                            mse[:,0:kmsetop+1,:,:], dp, g)

    # Horizontal advection
    hadv_dse_vint = hadv_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                            dse[:,0:kmsetop+1,:,:], dp, g)
    hadv_mse_vint = hadv_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                            mse[:,0:kmsetop+1,:,:], dp, g)

    # Mass divergence
    dse_diverg_vint = diverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                                dse[:,0:kmsetop+1,:,:], dp, g)
    mse_diverg_vint = diverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                                mse[:,0:kmsetop+1,:,:], dp, g)

    # Horizontal flux divergence
    dse_fluxdiverg_vint = tot_fdiverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                                        dse[:,0:kmsetop+1,:,:], dp, g)
    mse_fluxdiverg_vint = tot_fdiverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
                                        mse[:,0:kmsetop+1,:,:], dp, g)

    ### Add variables to list ##############################################

    var_list=[]
    var_list.append(pres)
    var_list.append(rho) # 3D
    var_list.append(theta_e) # 3D
    var_list.append(pw)
    var_list.append(pw_sat)
    var_list.append(vmf)
    var_list.append(vmfu)
    var_list.append(vmfd)
    var_list.append(condh)
    var_list.append(dse) # 3D
    var_list.append(mse) # 3D
    var_list.append(mse_vint)
    var_list.append(vadv_dse_vint)
    var_list.append(vadv_mse_vint)
    var_list.append(hadv_dse_vint)
    var_list.append(hadv_mse_vint)
    var_list.append(dse_diverg_vint)
    var_list.append(mse_diverg_vint)
    var_list.append(dse_fluxdiverg_vint)
    var_list.append(mse_fluxdiverg_vint)

    # del rho
    # del theta_e
    # del dse
    # del mse

    ### Write out variables ##############################################

    var_names, descriptions, units, dims_set = var_ncdf_metadata(dims3d=(nt,nz,nx1,nx2))

    write_ncfile(datdir+filename_out, var_list, var_names, descriptions, units, dims_set)

    return


##### Main loops and calculations ######################################

# Main read loops for 3D (dependent) variables

print('Running storm: ',storm)

for ktest in range(ntest):
# for ktest in range(1):

    test_str=tests[ktest]

    print()
    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    # for imemb in range(2,nmem):

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        print(datdir)

        start = runtimer()

        main_routine(datdir)

        end = runtimer()
        time_elapsed = end - start
        print("Time elapsed for member: ", time_elapsed)
