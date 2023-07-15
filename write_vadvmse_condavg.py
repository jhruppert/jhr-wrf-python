#!/usr/bin/env python
# coding: utf-8

# Diagnose vertical MSE advection after conducting conditional averaging by
# precip classification (and other), yielding condition-time-height-
# averaged fields, then write them out.
# 
# Also contains function read_vadvmse_condh to read these variables.
# 
# James Ruppert  
# jruppert@ou.edu  
# 7/13/23

from netCDF4 import Dataset
import numpy as np
import os
from thermo_functions import density_moist, esat, mixr_from_e
from write_ncfile import write_ncfile
from precip_class import precip_class
from cfads_functions import mask_edges
from time import time as runtimer
import sys


#### Main settings

# Relative to NCRF start
it_start = -12
it_end   = 24
it_start = -1
it_end   = 1

storm = 'haiyan'
# storm = 'maria'

filename_out='vadv_mse.nc' # this is for ALL variables in the var_names list

# msetop = 100 # top for MSE integrals

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    # tests = ['ctl','ncrf36h']
    # ALL TESTS
    tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    # tests = ['ctl','ncrf48h']
    # ALL TESTS
    tests = ['ctl','ncrf48h','crfon72h']#'ncrf36h']
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

if storm == 'haiyan':
    it0_ctl = 36
elif storm == 'maria':
    it0_ctl = 48

######################################################################

#### NetCDF variable read functions

def tidy_up(var):
    var = np.squeeze(var)
    var = np.ma.masked_invalid(var, copy=False)
    var = mask_edges(var)
    # Can swap in masking by TC track here
    var = np.ma.filled(var, fill_value=np.nan) # Abandon masking, replace mask with NaNs
    return var

def var_read(datdir, varname, it0, it1):
    varfil_main = Dataset(datdir+varname+'.nc')
    # var = varfil_main.variables[varname][it0:it1+1,...]
    var = varfil_main.variables[varname][...]
    varfil_main.close()
    var = tidy_up(var)
    return var

def qvar_read(datdir, it0, it1):
    qfile = Dataset(datdir+'q_int.nc')
    # q_int = qfile.variables['q_int'][:,it0:it1+1,...]
    q_int = qfile.variables['q_int'][...]
    qfile.close()
    q_int = tidy_up(q_int)
    return q_int

def mse_read(datdir, it0, it1):
    mse_file = Dataset(datdir+'mse_diag.nc')
    # mse = mse_file.variables['mse'][it0:it1+1,...]
    # rho = mse_file.variables['rho'][it0:it1+1,...]
    # vmfu = mse_file.variables['vmfu'][it0:it1+1,...]
    # vmfd = mse_file.variables['vmfd'][it0:it1+1,...]
    mse = mse_file.variables['mse'][...]
    rho = mse_file.variables['rho'][...]
    vmfu = mse_file.variables['vmfu'][...]
    vmfd = mse_file.variables['vmfd'][...]
    mse_file.close()
    rho = tidy_up(rho)
    mse = tidy_up(mse)
    return mse, rho, vmfu, vmfd

##### MSE & DSE functions ################################

# All calculations assume:
#   - Lagrangian convention, i.e., partial()/partial(t) + HADV + VADV = ...
#   - dp is a constant
#   - < > is vertical integral

def vadv(w, rho, invar, dp):
    g = 9.81 # m/s2
    omeg = w * (-1)*g*rho # omega from hydrostatic approximation
    vadv = omeg * np.gradient(invar, (-dp), axis=2)
    return vadv

# def hadv(u, v, x1d, y1d, invar):
#     adv_x = u * np.gradient(invar, x1d, axis=3)
#     adv_y = v * np.gradient(invar, y1d, axis=2)
#     hadv = adv_x + adv_y
#     return hadv

# def diverg(u, v, x1d, y1d, invar):
#     # Divergence of phi = < phi del . V >
#     dudx = np.gradient(u, x1d, axis=3) # /s
#     dvdy = np.gradient(v, y1d, axis=2) # /s
#     diverg = invar*(dudx + dvdy)
#     return diverg

######################################################################

#### NetCDF variable metadata

def var_ncdf_metadata():
    
    var_names = [
        'mse',
        'w',
        'vadv_mse',
        # 'hadv_mse',
        # 'mse_diverg',
    ]
    descriptions = [
        'moist static energy, calculated as cpT + gz + L_v*q',
        'vertical motion',
        'VADV of MSE',
        # 'HADV of MSE',
        # 'MSE*(del.V)',
    ]
    units = [
        'J/kg',
        'm/s',
        'J/kg/s',
        # 'J/kg/s',
        # 'J/kg/s',
    ]
    dims = ('navg', 'nt', 'nz')
    dim_names = [
        dims,
        dims,
        dims,
        # dims,
        # dims,
    ]

    return var_names, descriptions, units, dim_names

##### Conditional averaging (and weighting) ######################################

def get_condavg_settings():
    
    condavg_label = [
        'all',      # All unmasked points
        # 'non-precip',
        'deep',
        'cong',
        'shall',
        'strat',
        'anvil',
        'deepcong', # deep + cong
        'stratanv', # strat + anv
        'allrain',  # deep + cong + strat + anv
        # 'upward',   # upward-motion-weighted
        # 'downward', # downward-motion-weighted
        ]

    condavg_title = [
        'All',
        # 'Non Precip',
        'Dc',
        'Cg',
        'Sh',
        'St',
        'An',
        'Dc+Cg',
        'St+An',  
        'Dp+Cg+St+An',
        # 'Upward',
        # 'Downward',
        ]

    return condavg_label, condavg_title

condavg_label, condavg_title = get_condavg_settings()
ncond = len(condavg_label)

def conditional_avg(strat, vmfu, vmfd, var_stack):
    # Code modified from time_series_condavg.ipynb

    condavg_label, condavg_title = get_condavg_settings()
    ncond = len(condavg_label)
    shape = var_stack.shape
    nvar = shape[0]
    nt = shape[1]
    nz = shape[2]

    var_avg = np.zeros((ncond, nvar, nt, nz))

    # Internal functions
    def weighted_avg(var_stack, weights):
        num   = np.nansum(var_stack * weights, axis=(3,4))
        denom = np.nansum(weights,             axis=(3,4))
        return num/denom

    # Masks *OUT* according to the below conditions
    # Averages conditionally based on selected strat indices

    # Extend strat along a new dimension representing all variables to avoid loops in this task
    strat_extend = np.repeat(strat[:, np.newaxis, ...], nz, axis=1)
    strat_extend = np.repeat(strat_extend[np.newaxis, ...], nvar, axis=0)

    kcond=0
    # all = simple average over whole domain
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(3,4))

    # kcond+=1
    # # deep
    # condition = (strat_extend == 1)
    # var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    # kcond+=1
    # # congestus
    # condition = (strat_extend == 2)
    # var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    # kcond+=1
    # # shallow
    # condition = (strat_extend == 3)
    # var_avg[kcond, ...] = np.nanmean(var_stack, axis=(2,3), where=condition)

    for istrat in range(1,6):
        kcond+=1
        condition = (strat_extend == istrat)
        var_avg[kcond, ...] = np.nanmean(var_stack, axis=(3,4), where=condition)

    kcond+=1
    # deep + cong
    condition = ((strat_extend == 1) | (strat_extend == 2))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(3,4), where=condition)

    kcond+=1
    # strat + anv
    condition = ((strat_extend == 4) | (strat_extend == 5))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(3,4), where=condition)

    kcond+=1
    # allrain: deep + cong + strat + anv
    condition = ((strat_extend == 1) | (strat_extend == 2) | (strat_extend == 4) | (strat_extend == 5))
    var_avg[kcond, ...] = np.nanmean(var_stack, axis=(3,4), where=condition)

    # Weighting function

    # kcond+=1
    # # upward-weighted
    # # vmfu_extend = np.repeat(vmfu[np.newaxis, :, np.newaxis, ...], nvar, axis=0)
    # vmfu_extend = np.repeat(vmfu[:, np.newaxis, ...], nz, axis=1)
    # vmfu_extend = np.repeat(vmfu_extend[np.newaxis, ...], nvar, axis=0)
    # var_avg[kcond, ...] = weighted_avg(var_stack, weights=vmfu_extend)

    # kcond+=1
    # # downward-weighted
    # # vmfd_extend = np.repeat(vmfd[np.newaxis, ...], nvar, axis=0)
    # vmfd_extend = np.repeat(vmfd[:, np.newaxis, ...], nz, axis=1)
    # vmfd_extend = np.repeat(vmfd_extend[np.newaxis, ...], nvar, axis=0)
    # var_avg[kcond, ...] = weighted_avg(var_stack, weights=vmfd_extend)

    return var_avg

########################################################################
########################################################################
##### Main loops and calculations ######################################
########################################################################
########################################################################

# Main read loops for 3D (dependent) variables

print('Running storm: ',storm)

ntest=len(tests)
for ktest in range(ntest):

    test_str=tests[ktest]

    print()
    print('Running test: ',test_str)

    # Data read time steps
    #   FOR CTL (e.g., baseline = 36):
    #       it0 = baseline - it_start
    #       it1 = baseline + it_end
    #   FOR NCRF (and other sens. tests): IT0 = 0
    #       it0 = 0
    #       it1 = 0 + it_end
    it0 = 0
    it1 = it_end
    if test_str == 'ctl':
        it0 = it0_ctl + it_start
        it1 += it0_ctl

    # Loop over ensemble members

    for imemb in range(nmem):

        start = runtimer()

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        print(datdir)

        # READ VARIABLES

        q_int = qvar_read(datdir, it0, it1) # mm
        strat = precip_class(q_int)

        # 0: non-precipitating
        # 1: deep convective
        # 2: congestus
        # 3: shallow
        # 4: stratiform
        # 5: anvil (weaker rainfall)

        mse, rho, vmfu, vmfd = mse_read(datdir, it0, it1)
        nt,nz,nx1,nx2 = mse.shape

        # Winds
        # u = var_read(datdir,'U', it0, it1) # m/s
        # v = var_read(datdir,'V', it0, it1)
        w = var_read(datdir,'W', it0, it1)

        ### Add variables to list ##############################################

        var_list=[]
        var_list.append(mse)
        var_list.append(rho)
        # var_list.append(u)
        # var_list.append(v)
        var_list.append(w)
        var_stack = np.stack(var_list, axis=0)
        del var_list
        del mse
        del rho
        del w

        ### Run conditional averaging ##############################################

        var_avg = conditional_avg(strat, vmfu, vmfd, var_stack)

        ### Calculate advection/divergence ##############################################

        # Vertical advection
        # vadv_mse = vadv(w, rho, mse, dp)
        vadv_mse = vadv(var_avg[:,2,...], var_avg[:,1,...], var_avg[:,0,...], dp)

        # Horizontal advection
        # hadv_mse = hadv(u, v, x1d, y1d, mse)

        # Mass divergence
        # mse_diverg = diverg(u, v, x1d, y1d, mse)
        # del u
        # del v

        # # Horizontal flux divergence
        # dse_fluxdiverg_vint = tot_fdiverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
        #                                        dse[:,0:kmsetop+1,:,:], dp, g)
        # mse_fluxdiverg_vint = tot_fdiverg_vint(u[:,0:kmsetop+1,:,:], v[:,0:kmsetop+1,:,:], x1d, y1d,
        #                                        mse[:,0:kmsetop+1,:,:], dp, g)

        ### Make new list ##############################################

        var_list=[]
        var_list.append(var_avg[:,0,...]) # MSE
        var_list.append(var_avg[:,2,...]) # w
        var_list.append(vadv_mse)
        del var_avg
        del vadv_mse

        ### Write out variables ##############################################

        var_names, descriptions, units, dim_names = var_ncdf_metadata()

        write_ncfile(datdir+filename_out, var_list, var_names, descriptions, units, dim_names)

        end = runtimer()
        time_elapsed = end - start
        print("Time elapsed for member: ", time_elapsed)