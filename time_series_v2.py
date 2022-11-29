#!/usr/bin/env python
# coding: utf-8

# ### Notebook to plot TC tracks
# 
# James Ruppert  
# jruppert@ou.edu  
# 11/27/22


from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker, cm
import cartopy
import subprocess
import sys


# #### Variable selection

storm = 'haiyan'
# storm = 'maria'

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
elif storm == 'maria':
    tests = ['ctl','ncrf48h']

# How many members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 1

# Starting member to read
memb0=1

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 8 # radius (deg) limit for masking around TC center

# #### Directories

figdir = "/home/jamesrup/figures/tc/ens/tracks/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"

nums=np.arange(memb0,nmem+memb0,1)
nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)


##### Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/'+itest+'/'
varfil_main = Dataset(datdir+'post/d02/T.nc')
nz = varfil_main.dimensions['level'].size
nx1 = varfil_main.dimensions['lat'].size
nx2 = varfil_main.dimensions['lon'].size#-xmin-1
pres = varfil_main.variables['pres'][:] # hPa
varfil_main.close()

process = subprocess.Popen(['ls '+datdir+'/wrfout_d02_*'],shell=True,
    stdout=subprocess.PIPE,universal_newlines=True)
output = process.stdout.readline()
wrffil = output.strip() #[3]
varfil_main = Dataset(wrffil)
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()
lon1d=lon[0,:]
lat1d=lat[:,0]


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

# Account for crossing of Date Line
if storm == 'haiyan':
    offset = 180
else:
    offset = 0
lon_offset = dateline_lon_shift(lon, reverse=0)
lon_offset_plt = lon + lon_offset
lon_offset_plt -= offset


# Create figure with all tracks

# ### Plotting routines ##############################################

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

# select plotting area
plt_area=[lon1d[0], lon1d[-1], lat1d[0], lat1d[-1]] # W,E,S,N


# ### Combined plot ##############################################

# create figure
fig = plt.figure(figsize=(20,10))
proj = cartopy.crs.PlateCarree(central_longitude=offset)
# box_proj = ccrs.PlateCarree(central_longitude=0)
ax = fig.add_subplot(111,projection=proj)
# ax.set_title(ptrack + '-hPa RVor, Ens Memb '+str(imemb+1), fontsize=20)
ax.set_title(storm.capitalize()+' ('+ptrack + '-hPa RVor, 10Memb)', fontsize=20)

ax.set_prop_cycle(color=[
    '#1f77b4', '#1f77b4', '#aec7e8', '#aec7e8', '#ff7f0e', '#ff7f0e', '#ffbb78', '#ffbb78', '#2ca02c', '#2ca02c', '#98df8a', '#98df8a',
    '#d62728', '#d62728', '#ff9896', '#ff9896', '#9467bd', '#9467bd', '#c5b0d5', '#c5b0d5'])#,
    # '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    # '#17becf', '#9edae5'])#,
    # marker=["s"]*20)


for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'
    track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
    print(track_file)

    # Read track
    ncfile = Dataset(track_file)
    clon = ncfile.variables['clon'][:] # deg
    clat = ncfile.variables['clat'][:] # deg
    ncfile.close()
    nt = clon.shape[0]

    clon_shift = clon
    if storm == 'haiyan':
        clon_offset = dateline_lon_shift(clon, reverse=0)
        clon_shift += clon_offset
    clon_shift -= offset

    # storm track
    plt.plot(clon_shift, clat, linewidth=2, label=nustr[imemb])
    skip=24
    itim=np.arange(0,nt,skip)
    plt.plot(clon_shift[itim], clat[itim], "s")

# add map features
ax.add_feature(cartopy.feature.LAND,facecolor="lightgray") #land color
# ax.add_feature(cartopy.feature.OCEAN) #ocean color
ax.add_feature(cartopy.feature.COASTLINE)
# ax.add_feature(cartopy.feature.STATES)
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

ax.set_extent(plt_area)

plt.legend(loc="upper right")

# plt.show()
# plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'_'+memb_all[imemb]+'.png',dpi=200, facecolor='white', \
plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'.png',dpi=200, facecolor='white', \
            bbox_inches='tight', pad_inches=0.2)
plt.close()


# ### Single member plots ##############################################

for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'
    track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
    print(track_file)

    # Read track
    ncfile = Dataset(track_file)
    clon = ncfile.variables['clon'][:] # deg
    clat = ncfile.variables['clat'][:] # deg
    ncfile.close()
    nt = clon.shape[0]

    clon_shift = clon
    if storm == 'haiyan':
        clon_offset = dateline_lon_shift(clon, reverse=0)
        clon_shift += clon_offset
    clon_shift -= offset

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111,projection=proj)
    ax.set_title(ptrack + '-hPa RVor, Ens Memb '+str(imemb+1), fontsize=20)

    # storm track
    plt.plot(clon_shift, clat, linewidth=2, label=nustr[imemb])#, color='k')
    skip=24
    itim=np.arange(0,nt,skip)
    plt.plot(clon_shift[itim], clat[itim], "s", color='r')

    # add map features
    ax.add_feature(cartopy.feature.LAND,facecolor="lightgray") #land color
    # ax.add_feature(cartopy.feature.OCEAN) #ocean color
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.STATES)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    ax.set_extent(plt_area)

    # plt.show()
    plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'_'+memb_all[imemb]+'.png',dpi=200, facecolor='white', \
                bbox_inches='tight', pad_inches=0.2)
    plt.close()
