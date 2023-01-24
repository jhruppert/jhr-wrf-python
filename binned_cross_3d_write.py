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

from netCDF4 import Dataset
import numpy as np
import subprocess
import sys
from thermo_functions import theta_equiv, density_moist
from mask_tc_track import mask_tc_track


# #### Main settings

# IVAR: Index variable (independent var)
ivar_select = 'vmf'
# options: the, vmf

# BINVAR: Variable to bin (dependent var)
#   i.e., BINVAR will be averaged as a function of IVAR
binvar_select = 'the'
# options: vmf, the

# Calculate anomaly as deviation from xy-mean
do_prm_xy = 0
if (binvar_select == 'the'):
    do_prm_xy = 1

# Mask out all points except [stratiform/nonrain/etc], or switch off
nstrat=5 # istrat = -1, 0, 1, 2, 3
         # -1-all, 0-non-raining, 1-conv, 2-strat, 3-conv+strat

# Number of sample time steps
nt=12
# nt=2
hr_tag = str(np.char.zfill(str(nt), 2))


# #### Additional settings and directories

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'

# Time selection
# hr_tag = str(np.char.zfill(str(nt), 2))

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

# Shift starting-read time step for CRFON comparison
t0_test=0
if 'crfon' in tests[1]:
    t0_test=24 # CRFON is restarted at t=24 in NCRF
    # memb0=5 # for CRFFON test

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 6 # radius (deg) limit for masking around TC center


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


def var_read_3d(datdir,varname,t0,t1):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return var
def var_read_2d(datdir,varname,t0,t1):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return var


# Strat/Conv index subset
def get_strattag(istrat):
    if istrat == -1:
        strattag='all'
    elif istrat == 0:
        strattag='nonrain'
    elif istrat == 1:
        strattag='conv'
    elif istrat == 2:
        strattag='strat'
    elif istrat == 3:
        strattag='stratconv'
    return strattag


# #### NetCDF variable write function

def write_nc(file_out,nt,nz,nbins,pres,bin_axis,var_binned,ivar_mean):

    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    nz_dim = ncfile.createDimension('nz', nz)
    bin_dim = ncfile.createDimension('nbins', nbins-1)

    levs = ncfile.createVariable('pres', np.float64, ('nz',))
    levs.units = 'hPa'
    levs.long_name = 'pressure'
    levs[:] = pres

    binsv = ncfile.createVariable('bins', np.float64, ('nbins',))
    binsv.units = 'K'
    binsv.long_name = 'bin_axis'
    binsv[:] = bin_axis

    vmf_binned = ncfile.createVariable('vmf', np.float64, ('nt','nz','nbins',))
    vmf_binned.units = 'kg/s/K'
    vmf_binned.long_name = 'vertical mass flux summed over x,y'
    vmf_binned[:,:,:] = var_binned

    th_e = ncfile.createVariable('th_e', np.float64, ('nt','nz',))
    th_e.units = 'K'
    th_e.long_name = 'equiv potential temperature averaged over x,y'
    th_e[:,:] = ivar_mean

    ncfile.close()


# #### Index aka Bin variable

# Variable settings

# Theta-e (equivalent potential temperature)
if ivar_select == 'th_e':
    fmin=315; fmax=365 # K
    step=1
    bins=np.arange(fmin,fmax+step,step)
    ivar_tag = 'isent'
# Vertical mass flux
elif ivar_select == 'vmf':
    bins=np.logspace(-3,1.1,num=20)
    # bins=np.logspace(-3.5,0.7,num=20)
    bins=np.concatenate((-1.*np.flip(bins),bins))
    ivar_tag='vmf'

nbins = np.size(bins)

# Create axis of bin center-points for plotting
bin_axis = (bins[np.arange(nbins-1)]+bins[np.arange(nbins-1)+1])/2


# #### Main loops and compositing

# Main read loops for 3D (dependent) variables

ntest=2

for ktest in range(ntest):

    test_str=tests[ktest]

    # This has been tested for corresponding time steps:
    #   t0=37,1 are the first divergent time steps in CTL,NCRF
    #   t0=25,1 are the first divergent time steps in NCRF,CRFON
    if test_str == 'ctl':
        if tests[1] == 'ncrf36h':
            t0=36
        elif tests[1] == 'ncrf48h':
            t0=48
    elif test_str == 'ncrf36h':
        t0=t0_test
    elif test_str == 'ncrf48h':
        t0=t0_test
    elif 'crfon' in test_str:
        t0=0

    t0+=1 # add one time step since NCRF(t=0) = CTL
    t1 = t0+nt

    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    
        print('Running imemb: ',memb_all[imemb])
    
        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        datdir3d = datdir #+'v2/'
        print(datdir)

        # Localize to TC track
        # NOTE: Using copied tracking from CTL for NCRF tests
        # track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'
        trackfil_ex=''
        if 'crf' in test_str:
            trackfil_ex='_ctlcopy'
        track_file = datdir+'../../track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        # Required variables

        # Stratiform index
        varname = 'strat'
        strat = var_read_2d(datdir,varname,t0,t1) # 0-non-raining, 1-conv, 2-strat, 3-other/anvil

        # Index AKA Bin variable ("ivar")

        varname='T'
        tmpk = var_read_3d(datdir3d,varname,t0,t1) # K
        varname = 'QVAPOR'
        qv = var_read_3d(datdir3d,varname,t0,t1) # kg/kg

        # Theta-e (equivalent potential temperature)
        if ivar_select == 'th_e':
            ivar = theta_equiv(tmpk,qv,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
        # Vertical mass flux
        elif ivar_select == 'vmf':
            # Density
            rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3
            varname='W'
            ww = var_read_3d(datdir3d,varname,t0,t1) # m/s
            ivar = ww * rho # kg / m2 /s


        # Three-dimensional dependent variables ("var")

        # Vertical mass flux
        if binvar_select == 'vmf':
            varname='W'
            w = var_read_3d(datdir3d,varname,t0,t1) # m/s
            # Subtract area-average W
            # w_mn = np.mean(w,axis=(2,3))
            # w_mn_copy = np.repeat(np.repeat(w_mn[:,:,np.newaxis,np.newaxis], nx1, axis=2), nx2, axis=3)
            # w -= w_mn_copy
            rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3
            var = rho * w
        # Theta-e (equivalent potential temperature)
        elif binvar_select == 'the':
            var = theta_equiv(tmpk,qv,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

        ### Process and save variable ##############################################

        # Calculate var' as anomaly from x-y-average, using large-scale (large-radius) var avg
        if do_prm_xy == 1:
            # radius_ls=3
            # var_ls = mask_tc_track(track_file, radius_ls, var, lon, lat, t0, t1)
            var_ls = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
            var_ls_avg = np.ma.mean(var_ls,axis=(0,2,3))
            var -= var_ls_avg[np.newaxis,:,np.newaxis,np.newaxis]

        # Localize to TC track
        # var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
        ivar = mask_tc_track(track_file, rmax, ivar, lon, lat, t0, t1)

        # Normalization factor: equal for all classes
        # ncell = np.ma.MaskedArray.count(ivar[0,0,:,:])

        for istrat in range(-1,nstrat-1):
        # for istrat in range(nstrat-2,nstrat-1):

            var_tmp = np.copy(var)
            ivar_tmp = np.copy(ivar)

            # Mask out based on strat/conv
            if (istrat != -1) & (istrat != 3):
                # var_tmp = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), var_tmp, copy=True)
                ivar_tmp = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), ivar_tmp, copy=True)
            elif istrat == 3:
                ivar_tmp = np.ma.masked_where(((np.repeat(strat,nz,axis=1) == 0) | (np.repeat(strat,nz,axis=1) == 3)),
                    ivar_tmp, copy=True)

            ivar_mean = np.mean(ivar_tmp, axis=(2,3))

            # Replace masked elements with zeros or NaNs
            # var_tmp  = np.ma.filled(var_tmp, fill_value=0)
            # ivar_tmp = np.ma.filled(ivar_tmp, fill_value=np.nan)

            var_binned=np.zeros((nt,nz,nbins-1))

            # Bin the variables from (x,y) --> (bin)
            for it in range(nt):
                for ik in range(nz):
                    for ibin in range(nbins-1):
                        indices = (np.logical_and((ivar_tmp[it,ik,:,:] >= bins[ibin]), (ivar_tmp[it,ik,:,:] < bins[ibin+1]))).nonzero()
                    # Total across cells
                        var_binned[it,ik,ibin] = np.sum(var_tmp[it,ik,indices[0],indices[1]], dtype=np.float64)
                    # Total divided by n-all-cells
                        # var_binned[it,ik,ibin] = np.sum(var_tmp[it,ik,indices[0],indices[1]], dtype=np.float64) / ncell
                    # Mean across ID'd cells
                        # var_binned[it,ik,ibin] = np.mean(var_tmp[it,ik,indices[0],indices[1]], dtype=np.float64)

            # Write out to netCDF file
            strattag = get_strattag(istrat)
            ex_tag='t0'+str(t0)
            binvar_tag = binvar_select
            if (do_prm_xy == 1): binvar_tag+='_xyp'
            file_out = datdir+ivar_tag+'_'+binvar_tag+'_'+strattag+'_'+hr_tag+'hr_'+ex_tag+'.nc'
            write_nc(file_out,nt,nz,nbins,pres,bin_axis,var_binned,ivar_mean)