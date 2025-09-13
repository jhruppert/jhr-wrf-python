# Python script to run and save vortex tracking from WRF TC output
# 
# James Ruppert  
# jruppert@ou.edu  
# September 2022

from distutils.log import info
from netCDF4 import Dataset
import numpy as np
import sys
import subprocess
from relvort import *
from object_track import *
from read_functions import *
import sys


# Function to account for crossing of the Intl Date Line
def dateline_lon_shift(lon_in, reverse):
    if reverse == 0:
        lon_offset = np.zeros(lon_in.shape)
        lon_offset[np.where(lon_in < 0)] += 360
    else:
        lon_offset = np.zeros(lon_in.shape)
        lon_offset[np.where(lon_in > 180)] -= 360
    # return lon_in + lon_offset
    return lon_offset

def write_track_nc(file_out, nt, track, clon_offset):
    ncfile = Dataset(file_out,mode='w')

    time_dim = ncfile.createDimension('time', nt) # unlimited axis (can be appended to).

    clat = ncfile.createVariable('clat', np.float64, ('time',))
    clat.units = 'degrees_north'
    clat.long_name = 'clat'
    clat[:] = track[1,:]

    clon = ncfile.createVariable('clon', np.float64, ('time',))
    clon.units = 'degrees_east'
    clon.long_name = 'clon'
    clon[:] = track[0,:] + clon_offset

    ncfile.close()
    return


# Options
# istorm  = 'haiyan'
istorm  = 'maria'
# imemb   = 'memb_01'
# itest   = 'ctl'
# itest   = 'ncrf36h'
# itest   = 'ncrf48h'
# itest   = 'crfon60h'

# for itest in ['ctl', 'ncrf36h', 'STRATANVIL_ON', 'STRATANVIL_OFF', 'STRAT_OFF']:
# for itest in ['STRATANVIL_ON', 'STRATANVIL_OFF', 'STRAT_OFF']:
for itest in ['ctl', 'ncrf48h']:

    # ptrack  = 850 # tracking pressure level
    ptrack  = 600 # tracking pressure level
    # var_tag = 'rvor'
    # var_tag = 'avor'
    # var_tag = 'avor_850-600'
    for var_tag in ['avor', 'avor_850-600']:

        # ------------------------------------------

        # Ens members
        nmem = 10 # number of ensemble members (1-5 have NCRF)
        # nmem=1
        memb0=1
        nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
        nustr = np.char.zfill(nums, 2)
        memb_all=np.char.add('memb_',nustr)


        # For initializing tracks in sensitivity tests
        if itest == 'ctl':
            i_senstest=False
        else:
            i_senstest=True

        # Sensitivity tests basis and time step from that basis
        if itest == 'ncrf36h' or itest == 'STRATANVIL_ON' or itest == 'STRATANVIL_OFF' or itest == 'STRAT_OFF':
            test_basis='ctl'
            it_basis=36
        elif itest == 'ncrf48h':
            test_basis='ctl'
            it_basis=48
        elif itest == 'crfon60h':
            test_basis='ncrf36h'
            it_basis=24
        elif itest == 'crfon72h':
            test_basis='ncrf48h'
            it_basis=24
        else:
            test_basis=''

        #top = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
        top = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
        datdir2 = 'post/d02/'

        # Get dimensions
        datdir = top+istorm+'/'+memb_all[0]+'/'+itest+'/'+datdir2
        nt, nz, nx1, nx2, pres = get_file_dims(datdir)
        # pres = np.arange(1000,25,-25)

        # Get WRF file list
        datdir = top+istorm+'/'+memb_all[0]+'/'+itest+'/'
        wrffiles, lat, lon = get_wrf_filelist(datdir)

        # Level selection
        ikread = np.where(pres == ptrack)[0][0]

        for imemb in range(nmem):
        # for imemb in range(1):

            main = top+istorm+'/'+memb_all[imemb]+'/'+itest+'/'
            datdir = main+'post/d02/'
            print("Running ",main)

            # Prepare variable to use for tracking
            if var_tag == 'avor':

                track_file_tag = var_tag+'_'+str(round(pres[ikread]))+'hPa'

                # Horizontal wind
                # ufil = Dataset(datdir+'U.nc') # this opens the netcdf file
                # u = ufil.variables['U'][:,ikread,:,:] # m/s
                # ufil.close()
                # vfil = Dataset(datdir+'V.nc') # this opens the netcdf file
                # v = vfil.variables['V'][:,ikread,:,:] # m/s
                # vfil.close()
                fil = Dataset(datdir+'AVOR_HiRes.nc') # this opens the netcdf file
                var = fil.variables['AVOR'][:,ikread,:,:] # 10**-5 /s
                fil.close()
                # Calculate vorticity
                # var=relvort(u,v,lat1d,lon1d)

            elif var_tag == 'avor_850-600':

                track_file_tag = var_tag

                ikread = np.where((pres <= 850) & (pres >=600))[0]
                fil = Dataset(datdir+'AVOR_HiRes.nc') # this opens the netcdf file
                avor = fil.variables['AVOR'][:,ikread,:,:] # 10**-5 /s
                fil.close()
                var = np.mean(avor, axis=1)

            nt=np.shape(var)[0]

            if (lon.min() < 0) and (lon.max() > 0):
                lon_offset = dateline_lon_shift(lon, reverse=0)
            else:
                lon_offset = 0

            # Set basis starting point for tracking for sensitivity tests
            if i_senstest:
                track_file = main+'../'+test_basis+'/track_'+track_file_tag+'.nc'
                ncfile = Dataset(track_file)
                clon = ncfile.variables['clon'][it_basis] # deg
                clat = ncfile.variables['clat'][it_basis] # deg
                ncfile.close()
                basis = [clon, clat]
                itmax = 0
            else:
                if istorm == 'maria' and itest == 'ctl':
                    basis=[-49.7, 12.2]
                    itmax=60
                else:
                    basis=None
                    itmax=None

            # Run tracking
            track, f_masked = object_track(var, lon + lon_offset, lat, basis=basis, itmax=itmax)

            clon=track[0,:]
            clon_offset = dateline_lon_shift(clon, reverse=1)
            # clat=track[1,:]

            # Write out to netCDF file
            file_out = main+'track_'+track_file_tag+'.nc'
            write_track_nc(file_out, nt, track, clon_offset)

        print(f"Done running {var_tag} for {itest}!")
