#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out vertically integrated moisture variables.
# 
# James Ruppert  
# jruppert@ou.edu  
# 5/20/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()


# #### Main settings

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['ncrf36h']
    # tests = ['crfon60h']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h','ncrf36h']
    # tests = [tests[1],'crfon72h']

# Members
nmem = 10 # number of ensemble members
# nmem = 2
enstag = str(nmem)
memb0=1 # Starting member to read

nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

datdir2 = 'post/d02/'

# Get pressure
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
file = Dataset(datdir+'U.nc') # this opens the netcdf file
pres = file.variables['pres'][:] # hPa
file.close()
delp=(pres[0]-pres[1])*1e2 # Pa


# #### NetCDF variable read functions

def var_read_3d(datdir, varname, delp):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,:,:,:]
    varfil_main.close()
    # Vertical integral (in pressure)
    var = np.sum(var, axis=1)*(delp/9.81) # density is absorbed into the unit conversion; units: [mm]
    return var


# #### NetCDF variable write function

def write_vars(datdir,q_int):

    file_out = datdir+'q_int.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    nq = np.shape(q_int)[0]
    nq, nt, nx1, nx2 = q_int.shape

    q_dim    = ncfile.createDimension('nq', nq)
    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    x1_dim   = ncfile.createDimension('nx1', nx1)
    x2_dim   = ncfile.createDimension('nx2', nx2)

    q_int_nc = ncfile.createVariable('q_int', np.single, ('nq','nt','nx1','nx2',))
    q_int_nc.units = 'mm'
    q_int_nc.long_name = 'vertically integrated hydrometeor variables, calculated as 1/g(int)dp'
    q_int_nc[:,:,:,:] = q_int[:,:,:,:]

    ncfile.close()


# #### Main loops and calculations

ntest=len(tests)

for ktest in range(ntest):

    test_str=tests[ktest]

    print('Running test: ',test_str)

    # Loop over ensemble members

    # for imemb in range(nmem):
    imemb = comm.rank

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print(datdir)

    # Required variables

    # read in mixing ratios
    q_list = ['QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP']
    nvar = len(q_list)
    q_int = []
    for ivar in range(nvar):
        q_int.append(var_read_3d(datdir, q_list[ivar], delp))
    q_int = np.stack(q_int, axis=0)

    ### Write out variables ##############################################

    write_vars(datdir, q_int)