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

def object_track(f, lon, lat, nx_sm):

    shape=np.shape(f)
    nt,nx,ny = shape

    if len(shape) != 3:
        print("Check input dimensions!")
        print("Shape: ",shape)
        sys.exit()

    #############################################

    # SETTINGS
    # nx_sm=9    # Following Chavas 2013 (XX smoothing run x30 times)
    n_repeat=3 # Run smoothing n-times

    # MAXIMUM RADIUS RANGE FROM TIME-SPACE MAX
    r_max=5 # Masked out beyond this radius at the neighboring time steps [degrees]
        # ^ also used to determine buffer from boundaries
    dxkm=3 # Grid spacing [km]

    # CONSTANTS
    # m2deg=1./(111e3)
    # THE BELOW ARE FOR TRACKING VORTEX SUBJECT TO C_MAX
        # dt=(time[1]-time[0])*86400.d # s
        # idt=1d/dt
        # rmax_track = c_max * dt * m2deg # maximum single-time-step displacement [degrees]

    # 3-dimensional lon/lat for weighting
    lon_tim = np.repeat(lon[np.newaxis,:,:], nt, axis=0)
    lat_tim = np.repeat(lat[np.newaxis,:,:], nt, axis=0)

    # OBJECT TRACKING SUBJECT TO A MAXIMUM ALLOWABLE TRANSLATION SPEED
    # m2deg=1./(111e3)
    # c_max=40. # A maximum allowable translation speed (should be set to something very liberal)
    # dt=3600.  # Data time step [seconds]
    # idt=1./dt
    # rmax_track = c_max * dt * m2deg ; maximum single-time-step displacement [degrees]

    #############################################

    # Smooth input variable in x,y
    f_smooth = ndimage.uniform_filter(f,size=(0,nx_sm,nx_sm),mode='nearest')
    for ido in range(n_repeat-1):
        f_smooth = ndimage.uniform_filter(f_smooth,size=(0,nx_sm,nx_sm),mode='nearest')

    # Smooth input variable in time
    # n_repeat=3
    nt_smooth=3
    # Smooth input variable
    # for ido in range(n_repeat):
    f_smooth = ndimage.uniform_filter(f_smooth,size=(nt_smooth,0,0),mode='nearest')

    # Retain only values ≥ 3 sigma
    f_sigma = f_smooth / np.std(f_smooth)
    f_masked = np.ma.array(f_sigma)
    f_masked = np.ma.masked_where(np.abs(f_masked) < 3, f_masked) #copy=True)

    # Assurance
    # print(np.min(f_sigma))
    # print(np.max(f_sigma))
    # print(np.min(np.abs(f_sigma)))
    # print()
    # print(np.min(f_masked))
    # print(np.max(f_masked))
    # print(np.min(np.abs(f_masked)))

    #############################################

    # Track maxima in time by finding the centroid

    lon_tim = np.repeat(lon[np.newaxis,:,:], nt, axis=0)
    lat_tim = np.repeat(lat[np.newaxis,:,:], nt, axis=0)

    clon = np.average(lon_tim,axis=(1,2),weights=f_masked)
    clat = np.average(lat_tim,axis=(1,2),weights=f_masked)

    return f_masked
