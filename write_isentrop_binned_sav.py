#!/usr/bin/env python
# coding: utf-8

# ### Script to write out TC output binned according to a 3D variable to ncdf files.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 4/23/22


# NOTE: Using copied tracking from CTL for NCRF tests

# from netCDF4 import Dataset
import numpy as np
# import subprocess
import sys
from thermo_functions import *
from write_ncfile import *
from read_functions import *
from precip_class import *
from memory_usage import *
from mpi4py import MPI

# from mask_tc_track import mask_tc_track

comm = MPI.COMM_WORLD

print(comm.rank)
print()

# #### Main settings

storm = 'haiyan'
# storm = 'maria'

main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'

# Tests to read and compare
if storm == 'haiyan':
    # tests = ['ctl']
    # tests = ['ctl','ncrf36h']
    # tests = ['crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    tests = ['ctl','ncrf36h']
    tests = ['ncrf36h']
    tests = ['ctl','ncrf48h','ncrf36h']
    # tests = [tests[1],'crfon72h']
    # tests = ['crfon72h']
ntest = len(tests)

pclass_name = ['all','noncloud','deepc','congest','shallowc','strat','anvil']
npclass = len(pclass_name)

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 1

# TC tracking
# ptrack='600' # tracking pressure level
# var_track = 'rvor' # variable
# rmax = 6 # radius (deg) limit for masking around TC center
# rmax = 3

# Number of sample time steps
# nt=12
# nt=2
# hr_tag = str(np.char.zfill(str(nt), 2))

################################################
# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

# Get dimensions
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
nt_ctl, nz, nx1, nx2, pres = get_file_dims(datdir)

# Get WRF file list
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
wrffiles, lat, lon = get_wrf_filelist(datdir)
lat = lat[:,0]
lon = lon[0,:]

################################################
#### NetCDF variable metadata

def var_regrid_metadata(nt,nz,nbins,nbinsm1):

    var_names = [
        'bins',
        'pres',
        'theta_e_mn',
        'h_diabatic',
        'lw',
        'lwc',
        'sw',
        'swc',
        'w',
        'rho',
    ]
    descriptions = [
        'equivalent potential temperature bins',
        'pressure',
        'mean equivalent potential temperature',
        'H_DIABATIC',
        'LW heat tendency',
        'LW clear-sky heat tendency',
        'SW heat tendency',
        'SW clear-sky heat tendency',
        'vertical motion',
        'density',
    ]
    units = [
        'K',
        'hPa',
        'K',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'm/s',
        'kg/m^3',
    ]
    dims_all = (nt,nz,nbinsm1)
    dim_names = ('nt','nz','nbinsm1')
    dims_set = [
        [('nbins',),(nbins,)],
        [('nz',),(nz,)],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set

############################################################################

# #### Index aka Bin variable settings

# Theta-e (equivalent potential temperature)
fmin=305; fmax=365 # K
nbins = 70
bins=np.linspace(fmin,fmax,num=nbins)

# #### Main loops and compositing

# for ktest in range(ntest):
for ktest in range(1):

    test_str=tests[ktest]

    print()
    print('Running test: ',test_str)

    # Loop over ensemble members

    # for imemb in range(nmem):
    for imemb in range(1):

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2

        # Localize to TC track
        # NOTE: Using copied tracking from CTL for NCRF tests
        # track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'
        # trackfil_ex=''
        # if 'crf' in test_str:
        #     trackfil_ex='_ctlcopy'
        # track_file = datdir+'../../track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        # Required variables

        # Get dimensions
        nt, nz, nx1, nx2, pres = get_file_dims(datdir)
        buffer = 80
        nx1-=buffer*2
        nx2-=buffer*2

        t0=0
        nt=2
        t1=nt

        # Create single multi-dimensional variable to accommodate masking
        # nvars=8
        # all_vars = np.ma.zeros((nvars,nt,nz,nx1,nx2))
        # all_vars = np.zeros((nvars,nt,nz,nx1,nx2))

        # SPLITTING VARIABLES ACROSS CPUs
        # ALL NODES REQUIRE: PCLASS, THETA_E

        # Precip classification
        q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
        pclass = precip_class(q_int)
        pclass_z = np.repeat(pclass[:,np.newaxis,:,:], nz, axis=1)
        del q_int
        del pclass

        varname='T'
        tmpk = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K
        varname = 'QVAPOR'
        qv = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
        theta_e = theta_equiv(tmpk,qv,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

        if comm.rank != 6:
            del qv, tmpk

        if comm.rank == 0:
            varname = 'H_DIABATIC'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 1:
            varname = 'RTHRATLW'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 2:
            varname = 'RTHRATLWC'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 3:
            varname = 'RTHRATSW'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 4:
            varname = 'RTHRATSWC'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 5:
            varname = 'W'
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K/s
        elif comm.rank == 6:
            invar = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3
            del qv, tmpk

        ### Process and save variable ##############################################

        # Calculate var' as anomaly from x-y-average, using large-scale (large-radius) var avg
        # if do_prm_xy == 1:
        #     # radius_ls=3
        #     # var_ls = mask_tc_track(track_file, radius_ls, var, lon, lat, t0, t1)
        #     var_ls = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
        #     var_ls_avg = np.ma.mean(var_ls,axis=(0,2,3))
        #     var -= var_ls_avg[np.newaxis,:,np.newaxis,np.newaxis]

        # Localize to TC track
        # var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
        # ivar = mask_tc_track(track_file, rmax, ivar, lon, lat, t0, t1)

        # Normalization factor: equal for all classes
        # ncell = np.ma.MaskedArray.count(ivar[0,0,:,:])

        # for ipclass in range(npclass):
        for ipclass in range(1):

            print("Running ipclass: ",pclass_name[ipclass])

            # Mask out based on precipitation class
            if (ipclass > 0):
                indices = (pclass_z != ipclass)
                theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
                invar_masked = np.ma.masked_where(indices, invar, copy=True)
            else:
                # Create simple memory references
                # Delete these at end of loop so that originals aren't overwritten
                theta_e_masked = theta_e
                invar_masked = invar

            # Replace masked elements with zeros or NaNs
            # var_tmp  = np.ma.filled(var_tmp, fill_value=0)
            # ivar_tmp = np.ma.filled(ivar_tmp, fill_value=np.nan)

            if comm.rank ==0:

                theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))

                var_list_write=[]
                var_list_write.append(bins)
                var_list_write.append(pres)
                var_list_write.append(theta_e_mean)

            ############# SET UP SHARED MEMORY DISTRIBUTION #############

            # Code grabbed from https://stackoverflow.com/questions/32485122/shared-memory-in-mpi4py
            # on rank 0, create the shared block
            # on rank 1 get a handle to it (known as a window in MPI speak)
            dims = (7,nt,nz,nbins-1)
            size = np.prod(dims)
            itemsize = MPI.DOUBLE.Get_size() 
            if comm.rank == 0:
                nbytes = size * itemsize
            else:
                nbytes = 0

            # on rank 0, create the shared block
            # on rank 1 get a handle to it (known as a window in MPI speak)
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 

            # create a numpy array whose data points to the shared mem
            buf, itemsize = win.Shared_query(0) 
            assert itemsize == MPI.DOUBLE.Get_size() 
            allvars_binned = np.ndarray(buffer=buf, dtype='d', shape=dims)
            allvars_binned[:] = np.nan

            #############################################################

            # Bin the variables from (x,y) --> (bin)

            nmin = 2 # minimum points to average

            for it in range(nt):
                for iz in range(nz):
                    for ibin in range(nbins-1):
                        indices = ((theta_e_masked[0,it,iz,:,:] >= bins[ibin]) & (theta_e_masked[0,it,iz,:,:] < bins[ibin+1])).nonzero()
                        # Mean across ID'd cells
                        if indices[0].shape[0] > nmin:
                            allvars_binned[comm.rank, it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

            del theta_e_masked, invar_masked

            # in process rank 1:
            # write the numbers 0.0,1.0,..,4.0 to the first 5 elements of the array
            # if comm.rank == 1:
            #     ary[:5] = np.arange(5)

            # wait in process rank 0 until process 1 has written to the array
            comm.Barrier() 

            # check that the array is actually shared and process 0 can see
            # the changes made in the array by process 1
            if comm.rank == 0: 
                print(allvars_binned[:,1,:,30])

            if comm.rank == 0:

                for inode in range(7):
                    var_list_write.append(allvars_binned[inode])

                # Write out to netCDF file

                strattag = pclass_name[ipclass]
                # ex_tag='t0'+str(t0)
                # avgvar_tag = avgvar_select
                # radstr = str(rmax)
                # if (do_prm_xy == 1): avgvar_tag+='_xyp'
                # file_out = datdir+'binned_'+main_tag+'_'+avgvar_tag+'_'+strattag+'_'+hr_tag+'hr_'+ex_tag+'_rmax'+radstr+'.nc'

                file_out = datdir+'binned_isentrop_'+strattag+'.nc'

                var_names, descriptions, units, dims_set = var_regrid_metadata(nt,nz,nbins,nbins-1)

                write_ncfile(file_out, var_list_write, var_names, descriptions, units, dims_set)