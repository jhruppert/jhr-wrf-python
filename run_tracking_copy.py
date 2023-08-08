# Python script to run and save vortex tracking from WRF TC output
# 
# James Ruppert  
# jruppert@ou.edu  
# September 2022

# USING THIS SCRIPT TO SIMPLY COPY CTL TO NCRF LOCATIONS

from distutils.log import info
from netCDF4 import Dataset
import numpy as np
import sys
import subprocess
from relvort import relvort
from object_track import object_track


# Choices
ptrack  = 600 # tracking pressure level
istorm  = 'haiyan'
# itest   = 'ctl'
# itest   = 'ncrf36h'
# itest   = 'crfon60h'
itest   = 'STRAT_OFF'
itest   = 'STRATANVIL_OFF'
itest   = 'STRATANVIL_ON'

# itest   = 'ctl'
# istorm  = 'maria'
# itest   = 'ncrf48h'
# itest   = 'crfon72h'

var_tag = 'rvor'

# ------------------------------------------

print('Tracking at:',ptrack,'hPa')

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
if itest == 'ncrf36h':
    test_basis='ctl'
    it_basis=36
elif (itest == 'ncrf48h'):
    test_basis='ctl'
    it_basis=48
elif (itest == 'STRAT_OFF') or (itest == 'STRATANVIL_ON') or (itest == 'STRATANVIL_OFF'):
    test_basis='ctl'
    it_basis=48
elif itest == 'crfon60h':
    # Haiyan
    test_basis='ctl'
    it_basis=24+36
elif itest == 'crfon72h':
    # Maria
    test_basis='ctl'
    it_basis=24+48
else:
    test_basis=''

top = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# LonLat
main = top+istorm+'/memb_01/ctl/'
process = subprocess.Popen(['ls '+main+'wrfout_d02_*'],shell=True,
    stdout=subprocess.PIPE,universal_newlines=True)
output = process.stdout.readline()
m1ctl = output.strip() #[3]
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
main = top+istorm+'/'+memb_all[0]+'/'+itest+'/'
datdir = main+'post/d02/'
fil = Dataset(datdir+'U.nc') # this opens the netcdf file
pres = fil.variables['pres'][:] # hPa
fil.close()

# Level selection
ikread = np.where(pres == ptrack)[0][0]


for imemb in range(nmem):
# for imemb in range(1):

    main = top+istorm+'/'+memb_all[imemb]+'/'+itest+'/'
    datdir = main+'post/d02/'
    print("Running ",main)

    # Prepare variable to use for tracking
    if var_tag == 'rvor':

    #     # Horizontal wind
        ufil = Dataset(datdir+'U10.nc') # this opens the netcdf file
        var = ufil.variables['U10'][:,0,:,:] # m/s
        ufil.close()
    #     # vfil = Dataset(datdir+'V.nc') # this opens the netcdf file
    #     # v = vfil.variables['V'][:,ikread,:,:] # m/s
    #     # vfil.close()
    #     fil = Dataset(datdir+'AVOR.nc') # this opens the netcdf file
    #     var = fil.variables['AVOR'][:,ikread,:,:] # 10**-5 /s
    #     fil.close()

    #     # Calculate vorticity
    #     # var=relvort(u,v,lat1d,lon1d)

    nt=np.shape(var)[0]

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

    if (lon.min() < 0) and (lon.max() > 0):
        lon_offset = dateline_lon_shift(lon, reverse=0)
    else:
        lon_offset = 0

    # Set basis starting point for tracking for sensitivity tests
    if i_senstest:
        track_file = main+'../'+test_basis+'/track_'+var_tag+'_'+str(round(pres[ikread]))+'hPa.nc'
        ncfile = Dataset(track_file,mode='r')
        clon_ctl = ncfile.variables['clon'][it_basis:nt+it_basis] # deg
        clat_ctl = ncfile.variables['clat'][it_basis:nt+it_basis] # deg
        ncfile.close()
        # basis = [clon, clat]

        # Write out to netCDF file
        # file_out = main+'track_'+var_tag+'_'+str(round(pres[ikread]))+'hPa.nc'
        file_out = main+'track_'+var_tag+'_ctlcopy_'+str(round(pres[ikread]))+'hPa.nc'
        ncfile = Dataset(file_out,mode='w')

        time_dim = ncfile.createDimension('time', nt) # unlimited axis (can be appended to).

        clat = ncfile.createVariable('clat', np.float64, ('time',))
        clat.units = 'degrees_north'
        clat.long_name = 'clat'
        clat[:] = clat_ctl

        clon = ncfile.createVariable('clon', np.float64, ('time',))
        clon.units = 'degrees_east'
        clon.long_name = 'clon'
        clon[:] = clon_ctl

        ncfile.close()

    else:
        basis=0

#     # Run tracking
#     track, f_masked = object_track(var, lon + lon_offset, lat, i_senstest, basis)
    
#     clon=track[0,:]
#     clon_offset = dateline_lon_shift(clon, reverse=1)
#     # clat=track[1,:]


#     # Write out to netCDF file
#     file_out = main+'track_'+var_tag+'_'+str(round(pres[ikread]))+'hPa.nc'
#     ncfile = Dataset(file_out,mode='w')

#     time_dim = ncfile.createDimension('time', nt) # unlimited axis (can be appended to).

#     clat = ncfile.createVariable('clat', np.float64, ('time',))
#     clat.units = 'degrees_north'
#     clat.long_name = 'clat'
#     clat[:] = track[1,:]

#     clon = ncfile.createVariable('clon', np.float64, ('time',))
#     clon.units = 'degrees_east'
#     clon.long_name = 'clon'
#     clon[:] = track[0,:] + clon_offset

#     ncfile.close()

# print("Done!")
