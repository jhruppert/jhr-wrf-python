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

def track_object(f, lon, lat):

    shape=np.shape(f)
    nt,nx,ny = shape[0:2]

    # SETTINGS
    # nx_sm=round(0.5*1./(dims.lat[1]-dims.lat[0])) # Half-degree smoothing in space
    nx_sm=9    # Following Chavas 2013 (smoothing run x30 times)
    nt_sm=3    # temporal smoothing (n time steps)
    lmin=1e3   # Minimum distance threshold [km] between storms

    # CONSTANTS
    m2deg=1./(111e3)
    # THE BELOW ARE FOR TRACKING VORTEX SUBJECT TO C_MAX
    dt=(time[1]-time[0])*86400.d # s
    idt=1d/dt
    rmax_track = c_max * dt * m2deg # maximum single-time-step displacement [degrees]

    # Smooth input variable
    return shape

