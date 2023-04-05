#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out density (moist)
#       from TC output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 4/4/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import subprocess
import sys
from thermo_functions import density_moist


# #### Main settings

storm = 'haiyan'
storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
    # tests = [tests[1],'crfon60h']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    # tests = [tests[1],'crfon72h']

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 2
enstag = str(nmem)
# Starting member to read
memb0=1


nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

datdir2 = 'post/d02/'

##### Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
datdir3d = datdir #+'v2/'
varfil_main = Dataset(datdir3d+'T.nc')
nz = varfil_main.dimensions['level'].size
# lat = varfil_main.variables['XLAT'][:][0] # deg
# lon = varfil_main.variables['XLONG'][:][0] # deg
nx1 = varfil_main.dimensions['lat'].size
nx2 = varfil_main.dimensions['lon'].size
pres = varfil_main.variables['pres'][:] # hPa
dp = (pres[0]-pres[1])*1e2 # Pa
varfil_main.close()

process = subprocess.Popen(['ls '+main+storm+'/'+memb_all[0]+'/'+tests[0]+'/wrfout_d02_*'],shell=True,
    stdout=subprocess.PIPE,universal_newlines=True)
output = process.stdout.readline()
wrffil = output.strip() #[3]
varfil_main = Dataset(wrffil)
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()

# ikread = np.where(pres == pres_top)[0][0]


# #### NetCDF variable read functions


def var_read(datdir,varname):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,:,:,:]
    varfil_main.close()
    return var


# #### NetCDF variable write function

def write_vars(datdir,nt,nz,nx1,nx2,rho):

    file_out = datdir+'density.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    rho_nc = ncfile.createVariable('rho', np.single, ('nt','nz','nx1','nx2',))
    rho_nc.units = 'J/kg'
    rho_nc.long_name = 'density , calculated as p/(RT-virtual)'
    rho_nc[:,:,:,:] = rho[:,np.newaxis,:,:]

    ncfile.close()


# #### Main loops and calculations

# Main read loops for 3D (dependent) variables

ntest=2

for ktest in range(ntest):

    test_str=tests[ktest]

    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    
        print('Running imemb: ',memb_all[imemb])
    
        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        print(datdir)

        # Required variables

        # Water vapor
        varname='QVAPOR'
        qv = var_read(datdir,varname) # kg/kg
        nt,nz,nx1,nx2 = qv.shape
        # Temperature
        varname='T'
        tmpk = var_read(datdir,varname) # K

        # Calculate density
        rho = density_moist(tmpk,qv,pres[np.newaxis,:,np.newaxis,np.newaxis,]*1e2)

        ### Write out variables ##############################################

        write_vars(datdir,nt,nz,nx1,nx2,rho)