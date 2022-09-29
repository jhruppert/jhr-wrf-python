#!/usr/bin/env python
# coding: utf-8
# 
# Calculate relative vorticity
# 
# Input:
#   u   = zonal wind (x,y) in m/s
#   v   = meridional ...
#   lon = longitude points (deg) as lon = lon(x)
#   lat = longitude points (deg) as lat = lat(y)
# 
# Returns: array (x,y) of relative vorticity (/s).
# 
# James Ruppert
# September 2022
# 

import numpy as np
import sys

def relvort(u, v, lon, lat, londim, latdim):
    
    a = 6371e3 # Earth radius, m
    inva = 1./a
    inv_cos = 1./np.cos(np.radians(lat))

    dudy = np.gradient(u,lat,axis=latdim) * inva
    dvdx = np.gradient(v,lon,axis=londim) * inv_cos[np.newaxis,:,np.newaxis] * inva

    # print("Shape of gradient variable:",np.shape(dvdx))
    vor = (dvdx - dudy)

    return vor




