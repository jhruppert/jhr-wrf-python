#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out saturation fraction (a vertically integrated variable)
#       from TC output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/8/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import subprocess
import sys
from thermo_functions import esat, mixr_from_e


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
dp = (pres[1]-pres[0])*1e2 # Pa
varfil_main.close()

process = subprocess.Popen(['ls '+main+storm+'/'+memb_all[0]+'/'+tests[0]+'/wrfout_d02_*'],shell=True,
    stdout=subprocess.PIPE,universal_newlines=True)
output = process.stdout.readline()
wrffil = output.strip() #[3]
varfil_main = Dataset(wrffil)
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()


# #### NetCDF variable read functions


def var_read(datdir,varname):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,:,:,:]
    varfil_main.close()
    return var


# #### NetCDF variable write function

def write_vars(datdir,nt,nx1,nx2,satfrac):

    file_out = datdir+'satfrac.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', 1)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    var_nc = ncfile.createVariable('vmfu', np.float64, ('nt','nz','nx1','nx2',))
    var_nc.units = '%'
    var_nc.long_name = 'saturation fraction, as r / rv (mixing ratio)'
    var_nc[:,:,:] = satfrac[:,np.newaxis,:,:]

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

        # Mixing ratio
        varname='QVAPOR'
        qv = var_read(datdir,varname) # kg/kg
        nt,nz,nx1,nx2 = qv.shape

        # Temp
        varname='T'
        tmpk = var_read(datdir,varname) # K

        # Calculate r-star (sat mixing ratio)
        e_sat = esat(tmpk) # Pa
        qv_sat = mixr_from_e(e_sat, (pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/kg
        print(np.max(qv_sat))
        sys.exit()

        # Vertically integrate
        g = 9.81
        dp = (pres[0]-pres[1])*1e2
        cons = dp/g
        iktop = np.where(pres == 100)[0][0]+1 # Integrate up to 100 hPa
        pw     = np.sum(qv[:,0:iktop,:,:], axis=1) * cons
        pw_sat = np.sum(qv_sat[:,0:iktop,:,:], axis=1) * cons

        satfrac = pw / pw_sat

        # Fill masked points with nan
        vmfu = np.ma.filled(vmfu, fill_value=np.nan)
        vmfd = np.ma.filled(vmfd, fill_value=np.nan)

        # Convert from units of heat tendency to rain rate
        lv0=2.5e6 # J/kg
        cp=1004.  # J/K/kg
        condh *= cp/lv0 # kg*K/m2/s --> kg/m2/s = mm/s (i.e., rain rate)
        # Rain rate is in mm/day so do the same for MP heating:
        condh *= 3600*24 # mm/s --> mm/d

        ### Write out variables ##############################################

        write_vars(datdir,nt,nx1,nx2,vmfu,vmfd,condh)