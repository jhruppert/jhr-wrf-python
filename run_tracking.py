# Python script to run and save vortex tracking from WRF TC output
# 
# James Ruppert  
# jruppert@ou.edu  
# September 2022

from netCDF4 import Dataset
import numpy as np
import sys
from relvort import relvort
from object_track import object_track


# Choices
ptrack  = 600 # tracking pressure level
istorm  = 'haiyan'
# imemb   = 'memb_01'
itest   = 'ctl'
var_tag = 'rvor'

# ------------------------------------------

print('Tracking at:',ptrack,'hPa')

# Ens members
nmem = 20 # number of ensemble members (1-5 have NCRF)
memb0=1
nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

for imemb in range(nmem):

    wrfenkf = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
    main = wrfenkf+istorm+'/'+memb_all[imemb]+'/'+itest+'/'
    datdir = main+'post/d02/'
    print("Running ",main)


    # LonLat
    m1ctl = wrfenkf+istorm+'/memb_01/ctl/wrfout_d02_2013-11-01_00:00:00'
    fil = Dataset(m1ctl) # this opens the netcdf file
    lon = fil.variables['XLONG'][:][0] # deg
    lon1d=lon[0,:]
    lat = fil.variables['XLAT'][:][0] # deg
    lat1d=lat[:,0]
    fil.close()
    llshape=np.shape(lon)
    nx = llshape[1]
    ny = llshape[0]

    # Pressure
    fil = Dataset(datdir+'U.nc') # this opens the netcdf file
    pres = fil.variables['pres'][:] # hPa
    fil.close()

    # Level selection
    ikread = np.where(pres == ptrack)[0][0]

    # Prepare variable to use for tracking
    if var_tag == 'rvor':

        # Horizontal wind
        ufil = Dataset(datdir+'U.nc') # this opens the netcdf file
        u = ufil.variables['U'][:,ikread,:,:] # m/s
        ufil.close()
        vfil = Dataset(datdir+'V.nc') # this opens the netcdf file
        v = vfil.variables['V'][:,ikread,:,:] # m/s
        vfil.close()

        # Calculate vorticity
        var=relvort(u,v,lat1d,lon1d)

    nt=np.shape(var)[0]

    # Run tracking
    track, f_masked = object_track(var, lon, lat)
    # clon=track[0,:]
    # clat=track[1,:]


    # Write out to netCDF file
    file_out = main+'track_'+var_tag+'_'+str(round(pres[ikread]))+'hPa.nc'
    ncfile = Dataset(file_out,mode='w')

    time_dim = ncfile.createDimension('time', nt) # unlimited axis (can be appended to).

    clat = ncfile.createVariable('clat', np.float64, ('time',))
    clat.units = 'degrees_north'
    clat.long_name = 'clat'
    clat[:] = track[1,:]

    clon = ncfile.createVariable('clon', np.float64, ('time',))
    clon.units = 'degrees_east'
    clon.long_name = 'clon'
    clon[:] = track[0,:]

    ncfile.close()
