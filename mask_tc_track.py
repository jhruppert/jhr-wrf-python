# Function to apply masking beyond a radius threshold WRT to a tracked
#   object, e.g., a TC.
# 
# INPUTS:
#       track_file - file path and name for netCDF file containing track as clon,clat
#           (generated by the program run_tracking.py)
#       rmax - radius limit [degrees lon/lat]
#       var - variable to mask as (time,z,x,y) (note, the order of x,y is arbitrary)
#       lon - array of longitude points as (x,y) [deg]
#       lat - " " latitude points
#       t0, t1 - bounding time indices of var, assuming shorter in time than clon/clat
# 
# RETURNS:
#       returns a masked array of identical shape to var.
# 
# James Ruppert  
# jruppert@ou.edu  
# September 2022

from netCDF4 import Dataset
import numpy as np

def mask_tc_track(track_file, rmax, var, lon, lat, t0, t1):

    # Input dimensions
    nt,nz,nx1,nx2 = var.shape

    # Read TC track
    ncfile = Dataset(track_file)
    clon = ncfile.variables['clon'][:] # deg
    clat = ncfile.variables['clat'][:] # deg
    ncfile.close()

    # Calculate radius from center as array(time,x,y)
    lon3d = np.repeat(lon[np.newaxis,:,:], nt, axis=0)
    lat3d = np.repeat(lat[np.newaxis,:,:], nt, axis=0)
    lon3d -= clon[t0:t1,np.newaxis,np.newaxis]
    lat3d -= clat[t0:t1,np.newaxis,np.newaxis]
    radius3d = np.sqrt( lon3d**2 + lat3d**2 )
    
    # Add vertical dimension to match shape of var
    radius4d = np.repeat(radius3d[:,np.newaxis,:,:], nz, axis=1)

    # Apply mask
    var_mask = np.ma.masked_where(radius4d > rmax, var, copy=True)

    #XX Mask out domain edges

    return var_mask