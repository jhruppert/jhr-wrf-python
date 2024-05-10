#!/usr/bin/env python
# coding: utf-8

# ### Script to write out TC output binned according to a 3D variable to ncdf files.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 4/23/22

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

## NOTE: Using copied tracking from CTL for NCRF tests
# from mask_tc_track import mask_tc_track

# Parallelization notes:
#   Using mpi4py to distribute the work of reading and processing
#   large-dimensional numpy arrays. The processed results are then
#   passed back to the rank-0 node, which does the netcdf write-out.

# Testing mode:
#   all loops shortened to a single iteration and nt = 3
# testing=True
testing=False

comm = MPI.COMM_WORLD
nproc = comm.Get_size()



# #### Main settings

proc_var_list = ['tmpk', 'qv', 'rho', 'H_DIABATIC', 'RTHRATLW', 'RTHRATLWC', 'RTHRATSW', 'RTHRATSWC', 'W']
nvars = len(proc_var_list)

# Check for required number of processors
if nproc != nvars:
    print("Check NPROC (-n option)! Should be ",nvars)
    print("Killing batch job")
    sys.exit()

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

# Kill all loops after single iteration for testing mode
if testing:
    nmem = 1
    ntest = 1
    npclass = 1

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

def var_regrid_metadata(nt,nz,nbins):

    nbinsm1 = nbins-1

    var_names = [
        'bins',
        'pres',
        'theta_e_mn',
        'frequency',
        'tmpk',
        'qv',
        'rho',
        'h_diabatic',
        'lw',
        'lwc',
        'sw',
        'swc',
        'w',
    ]
    descriptions = [
        'equivalent potential temperature bins',
        'pressure',
        'mean equivalent potential temperature',
        'frequency',
        'temperature',
        'water vapor mixing ratio',
        'density',
        'H_DIABATIC',
        'LW heat tendency',
        'LW clear-sky heat tendency',
        'SW heat tendency',
        'SW clear-sky heat tendency',
        'vertical motion',
    ]
    units = [
        'K',
        'hPa',
        'K',
        'n-cells',
        'K',
        'kg/kg',
        'kg/m^3',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'm/s',
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
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set

############################################################################
# Driver functions
############################################################################

def get_pclass(datdir, t0, t1):
    # Precip classification
    q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
    pclass = precip_class(q_int)
    pclass_z = np.repeat(pclass[:,np.newaxis,:,:], nz, axis=1)
    return pclass_z

################################

def read_all_vars(datdir, t0, t1, proc_var_list):

    varname = 'theta_e'
    theta_e = read_mse_diag(datdir,varname,t0,t1,mask=True,drop=True) # K

    pclass_z = get_pclass(datdir,t0,t1)

    # Distribute variable processing onto all ranks.
    # Rank[0] then receives all processed results and does write-out.

    if comm.rank == 0:
        varname = 'T'
        invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K
    elif comm.rank == 1:
        varname = 'QVAPOR'
        invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
    elif comm.rank == 2:
        varname = 'rho'
        invar = read_mse_diag(datdir,varname,t0,t1,mask=True,drop=True) # kg/m3
    else:
        invar = var_read_3d(datdir,proc_var_list[comm.rank],t0,t1,mask=True,drop=True) # K/s or m/s

    return theta_e, pclass_z, invar

################################

def run_binning(ipclass, bins, theta_e, invar, pclass_z):

    # Mask out based on precipitation class
    if (ipclass > 0):
        indices = (pclass_z != ipclass)
        theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
        invar_masked = np.ma.masked_where(indices, invar, copy=True)
    else:
        # Create memory references
        # Delete these at end of loop so that originals aren't overwritten
        theta_e_masked = theta_e
        invar_masked = invar

    theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))

    # Bin the variables from (x,y) --> (bin)

    dims = (nt,nz,nbins-1)
    invar_binned = np.full(dims, np.nan)
    freq_binned = np.ndarray(dims, dtype=np.float64)

    nmin = 3 # minimum points to average

    if testing:
        nt_loop = 3
    else:
        nt_loop = nt

    for it in range(nt_loop):
        for iz in range(nz):
            for ibin in range(nbins-1):
                indices = ((theta_e_masked[it,iz,:,:] >= bins[ibin]) & (theta_e_masked[it,iz,:,:] < bins[ibin+1])).nonzero()
                binfreq = indices[0].size
                freq_binned[it,iz,ibin] = np.array(binfreq, dtype=np.float64)
                # Take mean across ID'd cells
                if binfreq > nmin:
                    invar_binned[it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

    return freq_binned, invar_binned, theta_e_mean

################################

def pclass_loop_write_ncdf(datdir, bins, theta_e, invar, pclass_z):

    for ipclass in range(npclass):

        if comm.rank == 0:
            print()
            print("Running ipclass: ",pclass_name[ipclass])

        freq_binned, invar_binned, theta_e_mean = run_binning(ipclass,bins,theta_e,invar,pclass_z)

        # Consolidate rebinned data onto Rank0 and write netCDF file

        if comm.rank > 0:

            comm.Send(np.ascontiguousarray(invar_binned, dtype=np.float64), dest=0, tag=comm.rank)

        else:

            var_list_write=[]
            var_list_write.append(bins)
            var_list_write.append(pres)
            var_list_write.append(theta_e_mean)
            var_list_write.append(freq_binned)
            var_list_write.append(invar_binned)

            for irank in range(1,nvars):
                dims = (nt,nz,nbins-1)
                invar_binned = np.empty(dims)
                comm.Recv(invar_binned, source=irank, tag=irank)
                # check that the unique arrays are appearing on process 0
                # print()
                # print(invar_binned[1,:,30])
                var_list_write.append(invar_binned)

            # Write out to netCDF file

            pclass_tag = pclass_name[ipclass]
            file_out = datdir+'binned_isentrop_'+pclass_tag+'.nc'
            var_names, descriptions, units, dims_set = var_regrid_metadata(nt,nz,nbins)
            write_ncfile(file_out, var_list_write, var_names, descriptions, units, dims_set)

    return

############################################################################
# Top-level loop
############################################################################

# #### Index aka Bin variable settings

# Theta-e (equivalent potential temperature)
fmin=305; fmax=365 # K
nbins = 70
bins=np.linspace(fmin,fmax,num=nbins)

# #### Main loops and compositing

for ktest in range(ntest):
# for ktest in range(1,2):

    test_str=tests[ktest]

    if comm.rank == 0:
        print()
        print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    # for imemb in range(1):

        if comm.rank == 0:
            print()
            print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2

        # Get dimensions
        nt, nz, nx1, nx2, pres = get_file_dims(datdir)

        # Account for dropped edges
        buffer = 80
        nx1-=buffer*2
        nx2-=buffer*2

        t0=0
        # nt=3
        t1=nt

        theta_e, pclass_z, invar = read_all_vars(datdir,t0,t1,proc_var_list)

        # Loop over pclass, bin data, write NetCDF
        pclass_loop_write_ncdf(datdir,bins,theta_e,invar,pclass_z)
