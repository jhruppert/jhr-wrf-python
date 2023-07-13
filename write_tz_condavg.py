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
from thermo_functions import density_moist, esat, mixr_from_e
from write_ncfile import write_ncfile
from time import time as runtimer
import sys



#### Main settings

storm = 'haiyan'
# storm = 'maria'

filename_out='condavg_tz_diag.nc' # this is for ALL variables in the var_names list

# msetop = 100 # top for MSE integrals

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
    # tests = ['crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    tests = ['ctl','ncrf48h']#'ncrf36h']
    # tests = [tests[1],'crfon72h']
    # tests = ['crfon72h']

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
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
varfil_main = Dataset(datdir3d+'T.nc')
nz = varfil_main.dimensions['level'].size
nx1 = varfil_main.dimensions['lat'].size
nx2 = varfil_main.dimensions['lon'].size
pres = varfil_main.variables['pres'][:] # hPa
dp = (pres[0]-pres[1])*1e2 # Pa
varfil_main.close()

# kmsetop = np.where(pres == msetop)[0][0]

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

def var_ncdf_metadata():
    
    var_names = [
        'mse',
        'vadv_mse',
        'hadv_mse',
        'mse_diverg',
        'mse_fluxdiverg',
    ]
    descriptions = [
        'moist static energy, calculated as cpT + gz + L_v*q',
        'VADV of MSE (up to 100 hPa)',
        'HADV of MSE (up to 100 hPa)',
        'MSE*(del.V) (up to 100 hPa)',
        # 'del.(MSE*V) (up to 100 hPa)',
    ]
    units = [
        'J/kg',
        'J/kg/s',
        'J/kg/s',
        'J/kg/s',
        # 'J/kg/s',
    ]
    # dims2d = ('nt','nx1','nx2')
    dims3d = ('nt','nz')
    dim_names = [
        dims3d,
        dims3d,
        dims3d,
        dims3d,
        # dims3d,
    ]
    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dim_names

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

# def vert_int(invar, dp, g):
#     # Vertically integrate: 1/g * SUM(var*dp)
#     # Negative is absorbed by dp>0
#     var_sum = np.sum(invar, axis=1)*dp/g
#     return var_sum

def vadv_vint(w, rho, invar, dp, g):
    omeg = w * (-1)*g*rho # omega from hydrostatic approximation
    vadv = omeg * np.gradient(invar, (-dp), axis=1)
    # vadv_vint = vert_int(vadv, dp, g)
    return vadv

def hadv_vint(u, v, x1d, y1d, invar):
    adv_x = u * np.gradient(invar, x1d, axis=3)
    adv_y = v * np.gradient(invar, y1d, axis=2)
    hadv = adv_x + adv_y
    # hadv_vint = vert_int((adv_x + adv_y), dp, g)
    return hadv

def diverg_vint(u, v, x1d, y1d, invar):
    # Divergence of phi = < phi del . V >
    dudx = np.gradient(u, x1d, axis=3) # /s
    dvdy = np.gradient(v, y1d, axis=2) # /s
    diverg = invar*(dudx + dvdy)
    # diverg = vert_int((invar * div), dp, g)
    return diverg

# def tot_fdiverg_vint(u, v, x1d, y1d, invar, dp, g):
#     # Total flux divergence of phi:
#     #   flux_diverg = < del . phi*V >
#     #       = < d/dx (phi*u) + d/dy (phi*v) >
#     diverg_x = np.gradient((u * invar), x1d, axis=3)
#     diverg_y = np.gradient((v * invar), y1d, axis=2)
#     flux_diverg = diverg_x + diverg_y
#     # flux_diverg_vint = vert_int(flux_diverg, dp, g)
#     return flux_diverg

##### Conditional averaging (and weighting) ######################################

def get_condavg_settings():
    
    condavg_label = [
        'all',      # All unmasked points
        # 'non-precip',
        'deep',
        'cong',
        'shall',
        # 'strat',
        # 'anvil',
        'deepcong', # deep + cong
        'stratanv', # strat + anv
        'allrain',  # deep + cong + strat + anv
        'upward',   # upward-motion-weighted
        'downward', # downward-motion-weighted
        ]

    condavg_title = [
        'All',
        # 'Non Precip',
        'Dc',
        'Cg',
        'Sh',
        # 'St',
        # 'An',
        'Dc+Cg',
        'St+An',  
        'Dp+Cg+St+An',
        'Upward',
        'Downward',
        ]

    return condavg_label, condavg_title

condavg_label, condavg_title = get_condavg_settings()
ncond = len(condavg_label)


def conditional_avg(strat, var_stack):

    condavg_label, condavg_title = get_condavg_settings()
    ncond = len(condavg_label)
    nvar = var_stack.shape[0]
    nt = np.shape(strat)[0]

    vmfu = var_stack[0]
    vmfd = var_stack[1]

    # var_avg = np.ma.zeros((ncond, nvar, nt))
    var_avg = np.zeros((ncond, nvar, nt))

    # Internal functions
    def weighted_avg(var_stack, weights):
        num   = np.nansum(var_stack * weights, axis=(2,3))
        denom = np.nansum(weights,             axis=(2,3))
        return num/denom

    # Masks *OUT* according to the below conditions
    # Averages conditionally based on selected strat indices

    # Extend strat along a new dimension representing all variables to avoid loops in this task
    strat_extend = np.repeat(strat[np.newaxis, ...], nvar, axis=0)

    kcond=0
    # all = simple average over whole domain
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3))

    kcond+=1
    # deep
    condition = (strat_extend == 1)
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    kcond+=1
    # congestus
    condition = (strat_extend == 2)
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    kcond+=1
    # shallow
    condition = (strat_extend == 3)
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    kcond+=1
    # deep + cong
    condition = ((strat_extend == 1) | (strat_extend == 2))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    kcond+=1
    # strat + anv
    condition = ((strat_extend == 4) | (strat_extend == 5))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    kcond+=1
    # allrain: deep + cong + strat + anv
    condition = ((strat_extend == 1) | (strat_extend == 2) | (strat_extend == 4) | (strat_extend == 5))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    # Weighting function

    kcond+=1
    # upward-weighted
    vmfu_extend = np.repeat(vmfu[np.newaxis, ...], nvar, axis=0)
    var_avg[kcond, ...] = weighted_avg(var_stack, weights=vmfu_extend)

    kcond+=1
    # downward-weighted
    vmfd_extend = np.repeat(vmfd[np.newaxis, ...], nvar, axis=0)
    var_avg[kcond, ...] = weighted_avg(var_stack, weights=vmfd_extend)

    return var_avg

##### Main loops and calculations ######################################

# Read base-state height once per storm
varname='ZB'
z_b = var_read(datdir,varname,nz) # m

# Main read loops for 3D (dependent) variables

print('Running storm: ',storm)

ntest=len(tests)
for ktest in range(ntest):
#for ktest in range(1,2):

    test_str=tests[ktest]

    print()
    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
#    for imemb in range(0,3):

        start = runtimer()

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        print(datdir)

        # Required variables

        # Water vapor
        varname='QVAPOR'
        qv = var_read(datdir,varname,nz) # kg/kg
        nt,nz,nx1,nx2 = qv.shape
        # Temperature
        varname='T'
        tmpk = var_read(datdir,varname,nz) # K
        # Height
        varname='Z'
        z = var_read(datdir,varname,nz) # m2/s2
        z += z_b

        # Dry static energy (DSE)
        g=9.81 # m/s2
        cp=1004. # J/K/kg
        dse = cp*tmpk + g*z # J/kg
        del z

        # Density
        rho = density_moist(tmpk,qv,pres[np.newaxis,:,np.newaxis,np.newaxis,]*1e2) # kg/m3

        # Calculate PW and PW_SAT for saturation fraction
        e_sat = esat(tmpk) # Pa
        qv_sat = mixr_from_e(e_sat, (pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/kg
        del e_sat
        # Integrate over full column
        pw = vert_int(qv, dp, g) # mm
        pw_sat = vert_int(qv_sat, dp, g) # mm
        del qv_sat

        # Latent heat of vaporization
        cp=1004.  # J/K/kg
        cpl=4186. # J/k/kg
        cpv=1885. # J/K/kg
        lv0=2.5e6 # J/kg
        lv = lv0 - (cpl-cpv)*(tmpk-273.15)
        del tmpk

        # Moist static energy (MSE)
        g=9.81 # m/s2
        mse = dse + lv*qv # J/kg
        del lv
        del qv
        mse_vint = vert_int(mse[:,0:kmsetop+1,:,:], dp, g) # J/m^2

        # Winds
        w = var_read(datdir,'W',nz) # m/s
        u = var_read(datdir,'U',nz)
        v = var_read(datdir,'V',nz)

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

        ### Additional diagnostics to save ##############################################

        # Microphysics (MP) latent heating
        lh = var_read(datdir,'H_DIABATIC',nz) # K/s
        condh = vert_int(lh[:,0:kmsetop,:,:], dp, g) # kg/*K/m2/s
        del lh

        # Vertical mass flux
        # Mask for Up/Dn
        wu = np.ma.masked_where((w < 0), w, copy=True)
        wd = np.ma.masked_where((w > 0), w, copy=True)
        del w
        vmfu = vert_int(wu[:,0:kmsetop,:,:], dp, g) # kg/m/s
        vmfd = vert_int(wd[:,0:kmsetop,:,:], dp, g) # kg/m/s
        del wu
        del wd

        ### Add variables to list ##############################################

        var_list=[]
        var_list.append(rho)
        var_list.append(pw)
        var_list.append(pw_sat)
        var_list.append(vmfu)
        var_list.append(vmfd)
        var_list.append(condh)
        var_list.append(dse)
        var_list.append(mse)
        var_list.append(mse_vint)
        var_list.append(vadv_dse_vint)
        var_list.append(vadv_mse_vint)
        var_list.append(hadv_dse_vint)
        var_list.append(hadv_mse_vint)
        var_list.append(dse_diverg_vint)
        var_list.append(mse_diverg_vint)
        var_list.append(dse_fluxdiverg_vint)
        var_list.append(mse_fluxdiverg_vint)

        del rho
        del dse
        del mse

        ### Write out variables ##############################################

        var_names, descriptions, units, dim_names = var_ncdf_metadata()

        write_ncfile(datdir+filename_out, var_list, var_names, descriptions, units, dim_names)

        end = runtimer()
        time_elapsed = end - start
        print("Time elapsed for member: ", time_elapsed)
