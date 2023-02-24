#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out moist static energy and its
#       vertical integral from TC output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/24/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import subprocess
import sys
# from thermo_functions import theta_equiv, density_moist


# #### Main settings

pres_top = 100

storm = 'haiyan'
# storm = 'maria'

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
nmem = 1#10 # number of ensemble members (1-5 have NCRF)
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

ikread = np.where(pres == pres_top)[0][0]


# #### NetCDF variable read functions


def var_read(datdir,varname,ikread):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,0:ikread+1,:,:]
    varfil_main.close()
    return var


# #### NetCDF variable write function

def write_vars(datdir,nt,nz,nx1,nx2,mse,mse_int):

    file_out = datdir+'mse.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    mse_nc = ncfile.createVariable('mse', np.float64, ('nt','nz','nx1','nx2',))
    mse_nc.units = 'J/kg'
    mse_nc.long_name = 'moist static energy, calculated as cpT + gz + L_v*q'
    mse_nc[:,:,:] = mse[:,np.newaxis,:,:]

    msei_nc = ncfile.createVariable('mse_int', np.float64, ('nt','nx1','nx2',))
    msei_nc.units = 'J/m^2'
    msei_nc.long_name = 'integrated moist static energy, calculated as 1/g*integral(mse)dp up to 100 hPa'
    msei_nc[:,:,:] = mse_int[:,np.newaxis,:,:]

    ncfile.close()


# #### Main loops and calculations

# Read base-state height once per storm
varname='ZB'
z_b = var_read(datdir,varname,ikread) # m

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
        qv = var_read(datdir,varname,ikread) # kg/kg
        nt,nz,nx1,nx2 = qv.shape
        # Temperature
        varname='T'
        tmpk = var_read(datdir,varname,ikread) # K
        # Height
        varname='Z'
        z = var_read(datdir,varname,ikread) # m2/s2
        z += z_b

        # Calculate Moist Static Energy (MSE)

        # ;LATENT HEAT OF VAPORIZATION
        cp=1004.  # J/K/kg
        cpl=4186. # J/k/kg
        cpv=1885. # J/K/kg
        lv0=2.5e6 # J/kg
        lv = lv0 - (cpl-cpv)*(tmpk-273.15)

        g=9.81 # m/s2
        mse = cp*tmpk + g*z + lv*qv

        # Vertically integrate
        cons = dp/g
        mse_int = np.sum(mse, axis=1) * cons # J/m^2

        ### Write out variables ##############################################

        write_vars(datdir,nt,nz,nx1,nx2,mse,mse_int)