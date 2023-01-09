#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out time-pressure data from TC output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/8/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import subprocess
import sys
from thermo_functions import theta_virtual
from mask_tc_track import mask_tc_track


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

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
# rmax = 6 # radius (deg) limit for masking around TC center

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

def write_vars(datdir,nt,nz,pres,var,vartag,var_units,var_longname):

    file_out = datdir+'time_pres_'+vartag+'.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)

    pres_nc = ncfile.createVariable('pres', np.float64, ('nz',))
    pres_nc.units = 'hPa'
    pres_nc.long_name = 'pressure'
    pres_nc[:] = pres[:]

    var_nc = ncfile.createVariable(vartag, np.float64, ('nt','nz',))
    var_nc.units = var_units
    var_nc.long_name = var_longname
    var_nc[:,:] = var[:,:]

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
    
        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'
        print(datdir)

        # track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
        # Localize to TC track
        # NOTE: Using copied tracking from CTL for NCRF tests
        trackfil_ex=''
        if 'ncrf' in tests[ktest]:
            trackfil_ex='_ctlcopy'
        track_file = datdir+'track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        datdir = datdir+datdir2

        # Required variables

        # Mixing ratio
        varname='AVOR'
        avor = var_read(datdir,varname) # 10^-5 /s
        nt,nz,nx1,nx2 = avor.shape

        # Mixing ratio
        varname='QVAPOR'
        qv = var_read(datdir,varname) # kg/kg

        # Temperature
        varname='T'
        tmpk = var_read(datdir,varname) # K

        # New calculations
        thv = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

        t0=0
        t1=nt

        # Run masking and averaging

        # AVOR
        rmax_avor = 1
        avor = mask_tc_track(track_file, rmax_avor, avor, lon, lat, t0, t1)
        avor_mn = np.mean(avor, axis=(2,3))

        # THV'

        # MESO-ALPHA (3-deg radius)
        rmax_alpha = 3
        thv_alpha = mask_tc_track(track_file, rmax_alpha, thv, lon, lat, t0, t1)
        thv_alpha_mn = np.mean(thv_alpha, axis=(2,3))

        # MESO-BETA (1-deg radius)
        rmax_beta = 1
        thv_beta = mask_tc_track(track_file, rmax_beta, thv, lon, lat, t0, t1)
        thv_beta_mn = np.mean(thv_beta, axis=(2,3))
        thv_prime = thv_beta_mn - thv_alpha_mn

        # Replace mask with NaN
        avor_mn = np.ma.filled(avor_mn, fill_value=np.nan)
        thv_prime = np.ma.filled(thv_prime, fill_value=np.nan)


        ### Write out variables ##############################################

        vartag='thvp'
        var_units='K'
        var_longname='virtual potential temperature anomaly (1 minus 3 deg radius)'
        write_vars(datdir,nt,nz,pres,thv_prime,vartag,var_units,var_longname)

        vartag='avor'
        var_units='10^-5 /s'
        var_longname='absolute vorticity'
        write_vars(datdir,nt,nz,pres,avor_mn,vartag,var_units,var_longname)