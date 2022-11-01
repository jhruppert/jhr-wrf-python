#!/usr/bin/env python
# coding: utf-8
# 
# Calculate relative vorticity
# 
# Input:
#   u   = zonal wind as (t,y,x) in m/s
#   v   = meridional ...
#   lon = longitude points (deg) as lon = lon(x)
#   lat = longitude points (deg) as lat = lat(y)
# 
# Returns: array (t,y,x) of relative vorticity (/s).
# 
# James Ruppert
# September 2022
# 

import numpy as np
import sys

def relvort(u, v, lat, lon):
    
    a = 6371e3 # Earth radius, m
    deg2rad = np.pi/180
    deg2meters = a * deg2rad
    cosf = np.cos(np.radians(lat))

    dudy = np.gradient( u , lat*deg2meters , axis=1)
    dvdx = np.gradient( v , lon*deg2meters , axis=2) / cosf[np.newaxis,:,np.newaxis]

    # print("Shape of gradient variable:",np.shape(dvdx))
    vor = (dvdx - dudy)

    return vor




