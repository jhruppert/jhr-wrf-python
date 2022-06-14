# ; 
# ; Track a TC or precursor vortex using an object-based algorithm, following
# ;   Davis et al. (2006, MWR) and Rios-Berrios et al. (2018, JAS).
# ; 
# ; The sell of this approach is that any 2D variable (vorticity, MSLP, precip)
# ;   can be used as input.
# ; 
# ; Steps:
# ;   1) Large-scale smoothing XXX in time and space.
# ;   2) Top-5% of variable is retained.
# ;   3) Centroid is found via weighted integral approach.
# ;   4) Sanity check for maximum possible phase speed for continuity.
# ;
# ; Returns: numpy array[itrack,2] where itrack corresponds to (potentially)
# ;   multiple identified tracks and the second dimension is (lon,lat).
# ;
# ; James Ruppert
# ; June 2022
# ; 

# SETTINGS
# nx_sm=round(0.5*1./(dims.lat[1]-dims.lat[0])) # Half-degree smoothing in space
nx_sm=9    # Following Chavas 2013 (smoothing run x30 times)
nt_sm=3    # temporal smoothing (n time steps)
lmin=1e3   # Minimum distance threshold [km] between storms
it_zero=24 # set all time steps prior to this one to NaN
c_max=20.  # a very liberal maximum allowable translation speed for tracking
pmin_thresh1=-4. # lifetime-minimum pressure as threshold for identifying a cyclone
# pmin_thresh1=-3.
pmin_thresh2=-3. # after storm ID, sequential time steps must retain center of this strength
pmin_thresh2=pmin_thresh1#-2.
nt_min=24. # minimum time steps for a cyclone to be retained

# CONSTANTS
m2deg=1./(111e3)
# THE BELOW ARE FOR TRACKING VORTEX SUBJECT TO C_MAX
dt=(time[1]-time[0])*86400.d # s
idt=1d/dt
rmax_track = c_max * dt * m2deg # maximum single-time-step displacement [degrees]

