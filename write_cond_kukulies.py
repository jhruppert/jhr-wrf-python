#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out condensation rate direct from WRF output using
# ### Julia Kukulies's functions.
# 
# James Ruppert  
# jruppert@ou.edu  
# 9/13/25


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
from read_functions import *
from mpi4py import MPI
# import subprocess
import PE_Kukulies_mpfunctions as micro
import wrf
import xarray as xr

comm = MPI.COMM_WORLD
nproc = comm.Get_size()


# #### Main settings

storm = 'haiyan'
storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl']
    # tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    # tests = ['ctl','ncrf48h','ncrf36h']
    tests = ['ctl','ncrf48h']

# Members
nmem = 10 # number of ensemble members
# nmem = 1
enstag = str(nmem)
# Starting member to read
memb0=1

nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

datdir2 = 'post/d02/'


# #### NetCDF variable write function

# def write_vars(datdir,theta_e):
def write_var(datdir, invar, filename, shortname, longname, units):

    file_out = datdir+filename
    ncfile = Dataset(file_out,mode='w', clobber=True)

    nt, nx1, nx2 = invar.shape

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    var_nc = ncfile.createVariable(shortname, np.double, ('nt','nx1','nx2',))
    var_nc.units = units
    var_nc.long_name = longname
    var_nc[:,:,:] = invar[:,:,:]

    ncfile.close()

# Read vars and get theta-e

# Function to get MSE
def get_condh(wrf_file):

    ds = xr.open_dataset(wrf_file)
    wrf_ds = Dataset(wrf_file)

    # get necessary variables
    vertical_velocity = wrf.getvar(wrf_ds, 'wa')
    temp = wrf.getvar(wrf_ds, 'tk')
    pressure = wrf.getvar(wrf_ds, 'pres')
    qcloud = ds.QCLOUD.squeeze()
    # base_pressure = ds.PB.squeeze()

    # output of get_condensation_rate is in kg/kg/s
    condensation_rate_t= micro.get_condensation_rate(vertical_velocity, temp, pressure)#, base_pressure)

    #create a cloud mask, because equation is conditional for grid cells where qcloud > 0 
    condensation_cloud = condensation_rate_t.where(qcloud > 0, 0) 
    # moreover, we are only interested in positive values, as negative values are evaporation 
    condensation_masked = condensation_cloud.where(condensation_cloud > 0, 0 ).data
    # integrate over pressure levels to get kg/m2/s
    condensation_rate = micro.pressure_integration(condensation_masked, -pressure.data)

    ds.close()
    wrf_ds.close()

    return condensation_rate

def run_job(datdir):

    # Get WRF files
    wrf_files, lat, lon = get_wrf_filelist(datdir)
    nfiles = len(wrf_files)

    ### Write out variables ##############################################

    # get variables
    cond_integrated = []
    for ifile in range(nfiles):
    # for ifile in range(2):
        icond = get_condh(wrf_files[ifile])
        cond_integrated.append(icond)
    cond_integrated = np.array(cond_integrated)

    metadata = ['condh_kukulies.nc', # filename
                'condh_k',           # shortname
                'p-integrated condensation rate based on Kukulies method', # longname
                'kg/m^2/s']          # units

    write_var(datdir+datdir2, cond_integrated, metadata[0], metadata[1], metadata[2], metadata[3])

    return None


# #### Main loops and calculations

ntest=len(tests)

for ktest in range(ntest):
# for ktest in range(1):

    test_str=tests[ktest]

    print('Running test: ',test_str)

    # Loop over ensemble members

    # if comm.rank == 0:
    #     imemb_select=[0,1,8]
    # elif comm.rank == 1:
    #     imemb_select=[2,3,9]
    # elif comm.rank == 2:
    #     imemb_select=[4,5]
    # elif comm.rank == 3:
    #     imemb_select=[6,7]
    # imemb_select=range(nmem)

    # nmemb_select = len(imemb_select)

    # for imemb in imemb_select:

    # Use a single node per ensemble member
    # imemb = comm.rank
    for ii in range(2):
        imemb = comm.rank + ii*5

    print()
    print('Rank: ',comm.rank)
    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'
    print(datdir)

    run_job(datdir)

    # Fix file names
    # process = subprocess.Popen(['mv '+datdir+'thetav_hires.nc '+datdir+'theta_e_HiRes.nc'],shell=True,
    # stdout=subprocess.PIPE,universal_newlines=True)
    # lines = process.stdout.readlines()
    # print(lines)

