# 
# Track a TC or precursor vortex using an object-based algorithm, following
#   Davis et al. (2006, MWR) and Rios-Berrios et al. (2018, JAS).
# 
# The sell of this approach is that any 2D variable (vorticity, MSLP, precip)
#   can be used as input.
# 
# Steps:
#   1) Large-scale smoothing XXX in time and space.
#   2) Top-5% of variable is retained.
#   3) Centroid is found via weighted integral approach.
#  xxx 4) Sanity check for maximum possible phase speed for continuity.
# 
# Input:
#   f   = input variable assumed to be in form f = f(t,x,y)
#   lon = longitude points (deg) as lon = lon(x)
#   lat = longitude points (deg) as lat = lat(y)
#  
# Returns: numpy array[itrack,2] where itrack corresponds to (potentially)
#   multiple identified tracks and the second dimension is (lon,lat).
# 
# James Ruppert
# June 2022
# 

import numpy as np
from scipy import ndimage
import sys

def object_track(f, lon, lat):

    shape=np.shape(f)
    nt,nx,ny = shape

    if len(shape) != 3:
        print("Check input dimensions!")
        print("Shape: ",shape)
        sys.exit()

    #############################################

    # SETTINGS

    nx_sm=9      # Following Chavas 2013 (XX smoothing run x30 times)
    nx_repeat=15 # Run x-smoothing n-times
    nt_smooth=3  # Time-step smoothing
    # nt_repeat=3  # Run time-smoothing n-times
    r_max=5      # Masked out beyond this radius [degrees] at adjacent time steps
                 # Also used to determine buffer from boundaries

    # 3-dimensional lon/lat for weighting
    lon3d = np.repeat(lon[np.newaxis,:,:], nt, axis=0)
    lat3d = np.repeat(lat[np.newaxis,:,:], nt, axis=0)

    #############################################

    # SMOOTHING

    # Smooth input variable in x,y
    f_smooth = ndimage.uniform_filter(f,size=(0,nx_sm,nx_sm),mode='nearest')
    for ido in range(nx_repeat-1):
        f_smooth = ndimage.uniform_filter(f_smooth,size=(0,nx_sm,nx_sm),mode='nearest')

    # Smooth input variable in time
    # for ido in range(nt_repeat):
    f_smooth = ndimage.uniform_filter(f_smooth,size=(nt_smooth,0,0),mode='nearest')

    # BULK MASKING

    # Mask out values < 3 sigma
    f_sigma = f_smooth / np.std(f_smooth)
    f_masked = np.ma.array(f_sigma)
    f_masked = np.ma.masked_where(np.abs(f_masked) < 3, f_masked, copy=False)

    # Mask out first time step
    f_masked[0,:,:] = np.nan
    f_masked = np.ma.masked_invalid(f_masked, copy=False)

    # Mask out data within 0.5*r_max from boundaries
    f_masked = np.ma.masked_where(lon3d <= lon1d[0]+0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lon3d >= lon1d[nx-1]-0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lat3d <= lat1d[0]+0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lat3d >= lat1d[ny-1]-0.5*r_max, f_masked, copy=False)

    #############################################

    # MASK BEYOND SPECIFIED RADIUS STARTING FROM ALL-TIME MAX

    # Locate the all-time maximum value
    fmax = np.max(f_masked)
    mloc=np.where(f_masked == fmax)
    itmax = mloc[0][0]
    xmax=mloc[2][0]
    ymax=mloc[1][0]

    radius = np.sqrt( (lon-lon1d[xmax])**2 + (lat-lat1d[ymax])**2 )

    # Mask surrounding points
    for it in range(itmax-1,np.minimum(itmax+1,nt-1)+1):
        f_masked[it,:,:] = np.ma.masked_where(radius > r_max, f_masked[it,:,:], copy=False)

    # Iterate downward from itmax
    for it in range(itmax-1,1,-1):

        fmax = np.max(f_masked[it,:,:])
        mloc = np.where(f_masked[it,:,:] == fmax)
        xmax = mloc[1][0]
        ymax = mloc[0][0]

        radius = np.sqrt( (lon-lon1d[xmax])**2 + (lat-lat1d[ymax])**2 )
        f_masked[it-1,:,:] = np.ma.masked_where(radius > r_max, f_masked[it-1,:,:], copy=False)

    # Iterate upward from itmax
    for it in range(itmax+1,nt-1):

        fmax = np.max(f_masked[it,:,:])
        mloc = np.where(f_masked[it,:,:] == fmax)
        xmax = mloc[1][0]
        ymax = mloc[0][0]

        radius = np.sqrt( (lon-lon1d[xmax])**2 + (lat-lat1d[ymax])**2 )
        f_masked[it+1,:,:] = np.ma.masked_where(radius > r_max, f_masked[it+1,:,:], copy=False)

    #############################################

    # TRACKING

    # Track maxima in time as the centroid
    clon = np.average(lon3d,axis=(1,2),weights=f_masked)
    clat = np.average(lat3d,axis=(1,2),weights=f_masked)

    track=np.stack([clon,clat])

    return track
