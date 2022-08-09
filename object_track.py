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

    # SETTINGS
    # nx_sm=9    # Following Chavas 2013 (XX smoothing run x30 times)
    n_repeat=3 # Run smoothing n-times
    # lmin=1e3   # Minimum distance threshold [km] between storms

    # CONSTANTS
    m2deg=1./(111e3)
    # THE BELOW ARE FOR TRACKING VORTEX SUBJECT TO C_MAX
        # dt=(time[1]-time[0])*86400.d # s
        # idt=1d/dt
        # rmax_track = c_max * dt * m2deg # maximum single-time-step displacement [degrees]

    # Smooth input variable
    f_smooth = ndimage.uniform_filter(f,size=(0,nx_sm,nx_sm),mode='nearest')
    for ido in range(n_repeat-1):
        f_smooth = ndimage.uniform_filter(f_smooth,size=(0,nx_sm,nx_sm),mode='nearest')

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

    return f_masked
