#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out precip classification NetCDF files.
# 
# James Ruppert  
# jruppert@ou.edu  
# 12/1/25


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
from mpi4py import MPI
from read_functions import *
from precip_class import *

# comm = MPI.COMM_WORLD
# nproc = comm.Get_size()


# Get dimensions
def get_nt(test_str):
    datdir = main+storm+'/'+memb_all[0]+'/'+test_str+'/'+datdir2
    nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
    # dp = (pres[1]-pres[0])*1e2 # Pa
    # nt=np.min([nt,nt_data-time_neglect])
    return nt_data, nz, pres

#### NetCDF variable write function
def write_vars(datdir,invar):

    file_out = datdir+'pclass.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    nt, nx1, nx2 = invar.shape

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    x1_dim   = ncfile.createDimension('nx1', nx1)
    x2_dim   = ncfile.createDimension('nx2', nx2)

    invar_nc = ncfile.createVariable('pclass', np.int_, ('nt','nx1','nx2',))
    invar_nc.units = '-'
    invar_nc.long_name = 'precip/cloud classification'
    invar_nc[:,:,:] = invar[:,:,:]

    ncfile.close()


# #### Main settings

# storm = 'haiyan'
# storm = 'maria'
for storm in ['haiyan','maria']:

    # main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
    main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

    # Tests to read and compare
    if storm == 'haiyan':
        tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
        # tests = ['ncrf36h']
        # tests = ['crfon60h']
    elif storm == 'maria':
        # tests = ['ctl','ncrf36h']
        tests = ['ctl','ncrf48h']#,'ncrf36h']
        # tests = [tests[1],'crfon72h']

    # Just do CTL
    tests = ['ctl']

    # Members
    nmem = 10 # number of ensemble members
    # nmem = 2
    enstag = str(nmem)
    memb0=1 # Starting member to read

    nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
    nustr = np.char.zfill(nums, 2)
    memb_all=np.char.add('memb_',nustr)

    datdir2 = 'post/d02/'


    # #### Main loops and calculations

    ntest=len(tests)

    for ktest in range(ntest):

        test_str=tests[ktest]

        print('Running test: ',test_str)

        # Loop over ensemble members

        # iimemb = comm.rank
        for imemb in range(nmem):
        # for ixmemb in range(2):

            # imemb = iimemb + 5*ixmemb
            # imemb = iimemb

            print('Running imemb: ',memb_all[imemb])

            datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
            print(datdir)
            print()

            # Required variables
            nt_test, nz, pres = get_nt(test_str)
            t1=nt_test
            q_int = read_qcloud(datdir,0,t1,mask=False,drop=False) # mm

            # Precip classification
            pclass = np.squeeze(precip_class(q_int))
            # pclass = pclass[:,np.newaxis,...] # Add single-element vertical dimension

            ### Write out variables ##############################################

            write_vars(datdir, pclass)
