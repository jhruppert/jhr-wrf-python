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
from read_functions import *
from thermo_functions import *
from mpi4py import MPI
# import subprocess


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
def write_var(datdir, invar, pres, filename, shortname, longname, units):

    file_out = datdir+filename
    ncfile = Dataset(file_out,mode='w', clobber=True)

    nt, nz, nx1, nx2 = invar.shape

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    var_pres = ncfile.createVariable('pres', np.single, ('nz',))
    var_pres.units = 'hPa'
    var_pres.long_name = 'pressure'
    var_pres[:] = pres[:]

    var_nc = ncfile.createVariable(shortname, np.single, ('nt','nz','nx1','nx2',))
    var_nc.units = units
    var_nc.long_name = longname
    var_nc[:,:,:,:] = invar[:,:,:,:]

    ncfile.close()

# Read vars and get theta-e

# Function to get MSE
def get_mse_dse(datdir, pres, t0, t1):
    # Constants
    cp=1004.  # J/K/kg
    cpl=4186. # J/k/kg
    cpv=1885. # J/K/kg
    lv0=2.5e6 # J/kg
    g=9.81 # m/s2

    # Required variables
    varname = 'T'
    tmpk = var_read_3d_hires(datdir,varname,t0,t1,mask=False) # K
    varname = 'QVAPOR'
    qv = var_read_3d_hires(datdir,varname,t0,t1,mask=False) # K
    z = var_read_3d_hires(datdir, 'Z', t0, t1, mask=False) # m

    # Read ZB once to add to Z
    zb = var_read_zb_hires(datdir, mask=False)
    z += zb

    # Latent heat of vaporization
    lv = lv0 - (cpl-cpv)*(tmpk-273.15)
    # Dry static energy (DSE)
    dse = cp*tmpk + g*z # J/kg
    # Moist static energy (MSE)
    mse = dse + lv*qv # J/kg

    return mse, dse

def run_job(datdir):

    # Get dimensions
    nt, nz, nx1, nx2, pres = get_file_dims(datdir)
    pres = np.arange(1000,25,-25)
    # nz = pres.shape[0]

    t0=0
    # nt=3
    t1=nt

    ### Write out variables ##############################################

    # get variables
    # theta_e, rho = get_theta_e_rho(datdir,pres,t0,t1)
    mse, dse = get_mse_dse(datdir,pres,t0,t1)

    metadata = ['mse_HiRes.nc', # filename
                'mse',          # shortname
                'moist static energy', # longname
                'J/kg']                # units

    write_var(datdir, mse, pres, metadata[0], metadata[1], metadata[2], metadata[3])

    metadata = ['dse_HiRes.nc', # filename
                'dse',          # shortname
                'dry static energy', # longname
                'J/kg']       # units

    write_var(datdir, dse, pres, metadata[0], metadata[1], metadata[2], metadata[3])

    return


# #### Main loops and calculations

ntest=len(tests)

for ktest in range(ntest):

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
    imemb = comm.rank

    print()
    print('Rank: ',comm.rank)
    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print(datdir)

    run_job(datdir)

    # Fix file names
    # process = subprocess.Popen(['mv '+datdir+'thetav_hires.nc '+datdir+'theta_e_HiRes.nc'],shell=True,
    # stdout=subprocess.PIPE,universal_newlines=True)
    # lines = process.stdout.readlines()
    # print(lines)

