# 
# Track a TC or precursor vortex using an object-based algorithm, following
#   Davis et al. (2006, MWR) and Rios-Berrios et al. (2018, JAS).
# 
# The sell of this approach is that any 2D variable (vorticity, MSLP, precip)
#   can be used as input.
# 
# Steps:
#   1) Large-scale smoothing XX in time and space.
#   2) Top-5% of variable is retained.
#   3) Centroid is found via weighted integral approach.
#  xxx 4) Sanity check for maximum possible phase speed for continuity.
# 
# Input:
#   f   = input variable assumed to be in form f = f(t,y,x)
#   lon = longitude points (deg) as lon = lon(y,x)
#   lat = longitude points (deg) as lat = lat(y,x)
#   sens_test = True or False, for if this case is a sensitivity test and
#       should have its tracking initiated on the basis of another simulation
#   basis = [lon, lat] the x,y location from which to begin tracking forward in
#       time, which is only used in the case of sens_test = True
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

def object_track(f, lon, lat, sens_test, basis):

    shape=np.shape(f)
    nt,ny,nx = shape

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
    # f_sigma = f_smooth / np.std(f_smooth)
    # f_masked = np.ma.array(f_sigma)
    f_std = (f_smooth - np.nanmean(f_smooth)) / np.nanstd(f_smooth)
    f_masked = np.ma.array(f_std)
    f_masked = np.ma.masked_where(np.abs(f_masked) < 3, f_masked, copy=False)

    # Mask out data within 0.5*r_max from boundaries
    f_masked = np.ma.masked_where(lon3d <= lon[0,0]   +0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lon3d >= lon[0,nx-1]-0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lat3d <= lat[0,0]   +0.5*r_max, f_masked, copy=False)
    f_masked = np.ma.masked_where(lat3d >= lat[ny-1,0]-0.5*r_max, f_masked, copy=False)

    #############################################

    if sens_test:

        # Assuming a sensitivity test restarted from another simulation
        print('Assuming sensitivity test; using basis to initialize track')

        # Therefore, using that restart time step as the basis from which
        # to track forward in time.

        radius = np.sqrt( (lon-basis[0])**2 + (lat-basis[1])**2 )
        itmax = 0

    else:

        # Assuming a cold start; will identify all-time-max object magnitude
        # and track forward and backward from there.
        print('Assuming cold start and using no basis')

        # Mask out first time step
        f_masked.mask[0,:,:] = True

        # Locate the all-time maximum value
        fmax = np.max(f_masked)
        mloc=np.where(f_masked == fmax)
        itmax = mloc[0][0]
        xmax=mloc[2][0]
        ymax=mloc[1][0]

        radius = np.sqrt( (lon-lon[ymax,xmax])**2 + (lat-lat[ymax,xmax])**2 )

    # Mask beyond specified radius at 2 neighboring time steps
    for it in range( np.maximum(itmax-1,0) , np.minimum(itmax+1,nt-1)+1 ):

        f_masked[it,:,:] = np.ma.masked_where(radius > r_max, f_masked[it,:,:], copy=False)

    # Do the same iterating all the way backward from itmax
    for it in range( np.maximum(itmax-1,0), 0, -1):

        fmax = np.max(f_masked[it,:,:])
        mloc = np.where(f_masked[it,:,:] == fmax)
        xmax = mloc[1][0]
        ymax = mloc[0][0]

        radius = np.sqrt( (lon-lon[ymax,xmax])**2 + (lat-lat[ymax,xmax])**2 )
        f_masked[it-1,:,:] = np.ma.masked_where(radius > r_max, f_masked[it-1,:,:], copy=False)

    # Do the same iterating all the way forward from itmax
    for it in range(itmax+1,nt-1):

        fmax = np.max(f_masked[it,:,:])
        mloc = np.where(f_masked[it,:,:] == fmax)
        xmax = mloc[1][0]
        ymax = mloc[0][0]

        radius = np.sqrt( (lon-lon[ymax,xmax])**2 + (lat-lat[ymax,xmax])**2 )
        f_masked[it+1,:,:] = np.ma.masked_where(radius > r_max, f_masked[it+1,:,:], copy=False)

    #############################################

    # TRACKING

    # Track maxima in time as the centroid = mean of latitude/longitude weighted by f
    clon = np.average(lon3d,axis=(1,2),weights=f_masked)
    clat = np.average(lat3d,axis=(1,2),weights=f_masked)

    track = np.ma.concatenate([clon[np.newaxis,:],clat[np.newaxis,:]])

    return track, f_masked
