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

# Testing mode: all loops shortened
testing=True
testing=False

comm = MPI.COMM_WORLD
nproc = comm.Get_size()


# #### Main settings

# proc_var_list = ['tmpk', 'qv', 'rho', 'H_DIABATIC', 'RTHRATLW', 'RTHRATLWC', 'RTHRATSW', 'RTHRATSWC', 'W']
proc_var_list = ['tmpk', 'theta_v', 'qv', 'rho', 'RTHRATLW', 'RTHRATLWC', 'RTHRATSW', 'RTHRATSWC', 'W']
nvars = len(proc_var_list)

# How many time steps to write out?
nt_write=24

# Use high vertical resolution output?
do_hires=True
# do_hires=False

storm = 'haiyan'
# storm = 'maria'

main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
postproc = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/jruppert/tc_postproc/"
datdir2 = 'post/d02/'

# Tests to read and compare
if storm == 'haiyan':
    # tests = ['ctl']
    # tests = ['ctl','ncrf36h']
    # tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    tests = ['ctl','ncrf36h']
    tests = ['ncrf36h']
    # tests = ['ctl','ncrf48h','ncrf36h']
    tests = ['ctl','ncrf48h']
    # tests = [tests[1],'crfon72h']
    # tests = ['crfon72h']
ntest = len(tests)

# pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','mcs','all']
pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','all']
npclass = len(pclass_name)

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)

# Kill all loops after single iteration for testing mode
if testing:
    nvars=1
    ntest = 1
    # nmem = 1
    nt_write=2
    # npclass = 1

################################################
# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

# Get dimensions
# datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
# nt_ctl, nz, nx1, nx2, pres = get_file_dims(datdir)

# # Get WRF file list
# datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
# wrffiles, lat, lon = get_wrf_filelist(datdir)
# lat = lat[:,0]
# lon = lon[0,:]

############################################################################
# Functions
############################################################################

#### NetCDF variable metadata

def var_regrid_metadata(nt,nz,nbins):

    nbinsm1 = nbins-1

    var_names = [
        'bins',
        'pres',
        'theta_e_mn',
        'pclass_frequency',
        'frequency',
        'tmpk',
        'tmpk_mean',
        'theta_v_prm',
        'thv_mean',
        'qv',
        'qv_mean',
        'rho',
        'rho_mean',
        # 'h_diabatic',
        'lw',
        'lw_mean',
        'lwc',
        'lwc_mean',
        'sw',
        'sw_mean',
        'swc',
        'swc_mean',
        'w',
        'w_mean',
    ]
    descriptions = [
        'equivalent potential temperature bins',
        'pressure',
        'mean equivalent potential temperature',
        'pclass frequency',
        'frequency',
        'temperature',
        'mean temperature',
        'virtual potential temperature xy-anomaly',
        'mean virtual potential temperature',
        'water vapor mixing ratio',
        'mean water vapor mixing ratio',
        'density',
        'mean density',
        # 'H_DIABATIC',
        'LW heat tendency',
        'mean LW heat tendency',
        'LW clear-sky heat tendency',
        'mean LW clear-sky heat tendency',
        'SW heat tendency',
        'mean SW heat tendency',
        'SW clear-sky heat tendency',
        'mean SW clear-sky heat tendency',
        'vertical motion',
        'mean vertical motion',
    ]
    units = [
        'K',
        'hPa',
        'K',
        'n-cells',
        'n-cells',
        'K',
        'K',
        'K',
        'K',
        'kg/kg',
        'kg/kg',
        'kg/m^3',
        'kg/m^3',
        # 'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'm/s',
        'm/s',
    ]
    dims_all = (nt,nz,nbinsm1)
    dim_names = ('nt','nz','nbinsm1')
    dims_set = [
        [('nbins',),(nbins,)],
        [('nz',),(nz,)],
        [('nt','nz'),(nt,nz)],
        [('nt',),(nt,)],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        # [dim_names,dims_all],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set

################################

def get_pclass(datdir, t0, t1):
    # Precip classification
    q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
    pclass = precip_class(q_int)
    pclass_z = np.repeat(pclass[:,np.newaxis,:,:], nz, axis=1)
    return pclass_z

################################

def get_theta_v(datdir,t0,t1,pres):
    varname = 'T'
    if do_hires:
        tmpk = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # K
    else:
        tmpk = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K
    varname = 'QVAPOR'
    if do_hires:
        qv = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
    else:
        qv = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
    theta_v = theta_virtual(tmpk, qv, pres[np.newaxis, :, np.newaxis, np.newaxis]*1e2)

    return theta_v

################################

def get_theta_v_prm(datdir,t0,t1,pres):

    theta_v = get_theta_v(datdir,t0,t1,pres)
    theta_v_mn = np.mean(theta_v, axis=(2,3))
    theta_v -= theta_v_mn[:, :, np.newaxis, np.newaxis]

    return theta_v

################################

def read_all_vars(datdir, t0, t1, ivar, proc_var_list):

    # Distribute variable processing onto all ranks.
    # Rank[0] then receives all processed results and does write-out.

    if ivar == 0:
        varname = 'T'
        if do_hires:
            invar = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # K
        else:
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # K
    elif ivar == 1:
        # Calculating theta-v is memory intensive since it requires
        # tmpk and qv, so do this first
        # invar = get_theta_v(datdir,t0,t1,pres)
        invar = get_theta_v_prm(datdir,t0,t1,pres)
    elif ivar == 2:
        varname = 'QVAPOR'
        if do_hires:
            invar = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
        else:
            invar = var_read_3d(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
    elif ivar == 3:
        varname = 'rho'
        if do_hires:
            # invar = density_moist(tmpk, qv, pres[np.newaxis, :, np.newaxis, np.newaxis]*1e2)
            invar = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # kg/kg
        else:
            invar = read_mse_diag(datdir,varname,t0,t1,mask=True,drop=True) # kg/m3
    else:
        if do_hires:
            invar = var_read_3d_hires(datdir,proc_var_list[ivar],t0,t1,mask=True,drop=True) # K/s or m/s
        else:
            invar = var_read_3d(datdir,proc_var_list[ivar],t0,t1,mask=True,drop=True) # K/s or m/s

    return invar

################################

def run_binning(ipclass, bins, theta_e, invar, pclass_z):

    shape = theta_e.shape
    nt = shape[0]
    nz = shape[1]
    nbins = bins.size

    # Mask out based on precipitation class
    # pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','mcs','all']
    if ipclass <= 5:
        indices = (pclass_z != ipclass)
        theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
        invar_masked = np.ma.masked_where(indices, invar, copy=True)
    elif ipclass == 6:
    #     # MCS = mask out NONCLOUD and SHALLOW
    #     indices = ((pclass_z != 0) & (pclass_z != 3))
    #     theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
    #     invar_masked = np.ma.masked_where(indices, invar, copy=True)
    # elif ipclass == 7:
        # Unmasked for "ALL" category
        theta_e_masked = theta_e
        invar_masked = invar

    # Frequency of cloud-type vs. time
    pclass_count = np.ndarray(nt, dtype=np.float64)
    for it in range(nt):
        pclass_count[it] = np.ma.count(theta_e_masked[it,2,:,:])

    theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))
    invar_mean = np.ma.mean(invar_masked, axis=(2,3))

    # Bin the variables from (x,y) --> (bin)

    dims = (nt,nz,nbins-1)
    invar_binned = np.full(dims, np.nan)
    freq_binned = np.ndarray(dims, dtype=np.float64)

    nmin = 3 # minimum points to average

    for it in range(nt):
        for iz in range(nz):
            for ibin in range(nbins-1):
                indices = ((theta_e_masked[it,iz,:,:] >= bins[ibin]) & (theta_e_masked[it,iz,:,:] < bins[ibin+1])).nonzero()
                binfreq = indices[0].size
                freq_binned[it,iz,ibin] = np.array(binfreq, dtype=np.float64)
                # Take mean across ID'd cells
                if binfreq > nmin:
                    invar_binned[it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

    return freq_binned, invar_binned, theta_e_mean, invar_mean, pclass_count

################################

def run_binning_multi(ipclass, bins, theta_e, invar, pclass_z):

    shape = theta_e.shape
    nt = shape[0]
    nz = shape[1]
    nbins = bins.size

    # Mask out based on precipitation class
    # pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','mcs','all']
    if ipclass <= 5:
        indices = (pclass_z != ipclass)
        theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
        invar_masked = np.ma.masked_where(indices, invar, copy=True)
    elif ipclass == 6:
    #     # MCS = mask out NONCLOUD and SHALLOW
    #     indices = ((pclass_z != 0) & (pclass_z != 3))
    #     theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
    #     invar_masked = np.ma.masked_where(indices, invar, copy=True)
    # elif ipclass == 7:
        # Unmasked for "ALL" category
        theta_e_masked = theta_e
        invar_masked = invar

    # Frequency of cloud-type vs. time
    pclass_count = np.ndarray(nt, dtype=np.float64)
    for it in range(nt):
        pclass_count[it] = np.ma.count(theta_e_masked[it,2,:,:])

    theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))
    invar_mean = np.ma.mean(invar_masked, axis=(2,3))

    # Bin the variables from (x,y) --> (bin)

    dims = (nvars,nt,nz,nbins-1)
    invar_binned = np.full(dims, np.nan)

    nmin = 3 # minimum points to average

    for it in range(nt):
        for iz in range(nz):
            for ibin in range(nbins-1):
                indices = ((theta_e_masked[it,iz,:,:] >= bins[ibin]) & (theta_e_masked[it,iz,:,:] < bins[ibin+1])).nonzero()
                binfreq = indices[0].size
                # Take mean across ID'd cells
                if binfreq > nmin:
                    for kvar in range(nvars):
                        invar_binned[it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

    return all_vars_binned, all_vars_mean

################################

def driver_loop_write_ncdf(datdir, outdir, bins, dims, t0, t1, proc_var_list):

    nt = dims[0]
    nz = dims[1]

    # Read theta-e
    varname = 'theta_e'
    if do_hires:
        # varname = 'T'
        # tmpk = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # K
        # varname = 'QVAPOR'
        # qv = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # K
        # theta_e = theta_equiv(tmpk, qv, qv, pres[np.newaxis, :, np.newaxis, np.newaxis]*1e2)
        theta_e = var_read_3d_hires(datdir,varname,t0,t1,mask=True,drop=True) # K
    else:
        theta_e = read_mse_diag(datdir,varname,t0,t1,mask=True,drop=True) # K

    # Get PCLASS
    pclass_z = get_pclass(datdir,t0,t1)


    for ipclass in range(npclass):
    # for ipclass in range(0,2):

        print()
        print("Running ipclass: ",pclass_name[ipclass])

        var_list_write=[]
        var_list_write.append(bins)
        var_list_write.append(pres)

        tmpk = read_all_vars(datdir,t0,t1,0,proc_var_list)
        freq_binned, tmpk_binned, theta_e_mean, tmpk_mean, pclass_count = run_binning(ipclass,bins,theta_e,tmpk,pclass_z)

        var_list_write.append(theta_e_mean)
        var_list_write.append(pclass_count)
        var_list_write.append(freq_binned)
        var_list_write.append(tmpk_binned)
        var_list_write.append(tmpk_mean)

        # Read variables
        all_vars = []
        for ivar in range(1,nvars):
            invar = read_all_vars(datdir,t0,t1,ivar,proc_var_list)
            all_vars.append(invar)

            # freq_binned, invar_binned, theta_e_mean, invar_mean, pclass_count = run_binning(ipclass,bins,theta_e,invar,pclass_z)

            # var_list_write.append(invar_binned)
            # var_list_write.append(invar_mean)

        all_vars_binned, all_vars_mean = run_binning_multi(ipclass,bins,theta_e,all_vars,pclass_z)

        # Write out to netCDF file

        pclass_tag = pclass_name[ipclass]
        if do_hires:
            file_out = outdir+'binned_isentrop_'+test_str+'_'+pclass_tag+'_HiRes.nc'
        else:
            file_out = outdir+'binned_isentrop_'+test_str+'_'+pclass_tag+'.nc'
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
# for ktest in range(3,ntest):
# for ktest in range(1):

    test_str=tests[ktest]

    # if comm.rank == 0:
    print()
    print('Running test: ',test_str)

    # Loop over ensemble members

    # for imemb in range(nmem):
    # for imemb in range(6,7):

    imemb = comm.rank

    # Skip some members
    # if ktest == 3 & imemb < 2:
    #     continue

    # if comm.rank == 0:
    print()
    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    outdir = postproc+storm+'/'+memb_all[imemb]+'/'

    # Get dimensions
    nt_test, nz, nx1, nx2, pres = get_file_dims(datdir)
    nt = np.min([nt_test,nt_write])

    if do_hires:
        pres = np.arange(1000,25,-25)
        nz = pres.shape[0]

    t0=0
    # nt=4
    t1=nt
    if test_str == 'ctl':
        t0+=36
        t1+=36

    dims = (t1-t0,nz)

    driver_loop_write_ncdf(datdir, outdir, bins, dims, t0, t1, proc_var_list)
