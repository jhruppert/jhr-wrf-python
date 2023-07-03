#!/usr/bin/env python
# coding: utf-8

# Python script to read in variables according to a list, smooth them, and
# write them out to a new file.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/24/23

from netCDF4 import Dataset
import numpy as np
import xarray as xr
from write_ncfile import write_ncfile
import os
import sys
from time import time as runtimer


#### Main settings

# Size of running-mean smoothing window (nx, ny)
dx=3 # km
xwindow=30 # km
nwindow=int(xwindow/dx)

storm = 'haiyan'
# storm = 'maria'

filename_read='mse_diag.nc' # this is for ALL variables in the var_names list
filename_out='msevars_smooth.nc' # this is for ALL MSE variables

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
    # tests = ['STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    # tests = [tests[1],'crfon72h']

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

# Get pressure
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
varfil_main = Dataset(datdir+'T.nc')
pres = varfil_main.variables['pres'][:] # hPa
varfil_main.close()

# WRFOUT file list to get Lat/Lon
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
dirlist = os.listdir(datdir)
subs="wrfout_d02"
wrf_files = list(filter(lambda x: subs in x, dirlist))
wrf_files.sort()
wrf_files = [datdir + s for s in wrf_files]
ncfile = Dataset(wrf_files[0])
lat = ncfile.variables['XLAT'][0,:,0]
lon = ncfile.variables['XLONG'][0,0,:]
ncfile.close()

# Address dateline issue
lon[(lon < 0)] += 360

##### Main loops and calculations ######################################

# var_names, long_names, units, dim_names = var_ncdf_metadata()

var_names = [
    'rho',
    'pw',
    'pw_sat',
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

# Main read loops for 3D (dependent) variables

ntest=len(tests)
for ktest in range(ntest):
# for ktest in range(0,1):

    test_str=tests[ktest]

    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):

        start = runtimer()

        print('Running imemb: ',memb_all[imemb])
    
        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2

        # Variable list to process

        var_smooth_list=[]
        descriptions=[]
        units=[]
        dim_names=[]

        for ivar in range(len(var_names)):
            ncfile = Dataset(datdir+filename_read)
            var = ncfile.variables[var_names[ivar]]
            descriptions.append(var.description)
            units.append(var.units)
            dim_names.append(var.dimensions)
            var = var[...]
            ncfile.close()

            # Variable settings
            nt=var.shape[0]
            time=np.arange(nt)

        ### Smooth variable ##############################################

            da = xr.DataArray(data=var, dims=dim_names[ivar])

            # print("STARTING SMOOTHING")
            # start = runtimer()
            var_smooth = da.rolling({dim_names[ivar][-2]:nwindow, dim_names[ivar][-1]:nwindow}, center=True).mean()
            # print("DONE SMOOTHING")
            # end = runtimer()
            # print(end - start) # Time in seconds, e.g. 5.38091952400282

        ### Add variable to list ##############################################

            var_smooth_list.append(var_smooth.data)

        ### Write out variables ##############################################

        write_ncfile(datdir+filename_out, var_smooth_list, var_names, descriptions, units, dim_names)

        end = runtimer()
        time_elapsed = end - start
        print("Time elapsed for member: ", time_elapsed)