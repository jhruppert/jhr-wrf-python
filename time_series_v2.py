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
from mask_tc_track import mask_tc_track
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

# Strat/Conv index subset
istrat=2

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 8 # radius (deg) limit for masking around TC center

# #### Directories

figdir = "/home/jamesrup/figures/tc/ens/tracks/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

nums=np.arange(memb0,nmem+memb0,1)
nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)


# Strat/Conv index subset
if istrat == -1:
    fig_extra=''
else:
    if istrat == 0:
        strattag='Nonrain'
    elif istrat == 1:
        strattag='Conv'
    elif istrat == 2:
        strattag='Strat'
    elif istrat == 3:
        strattag='Anv'
    fig_extra='_'+strattag.lower()

##### Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
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


# Create figure with all tracks

# ### Plotting routines ##############################################

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


# ### Combined plot ##############################################

# create figure
fig = plt.figure(figsize=(14,5))
ax = fig.add_subplot(111)

ax.set_title(storm.capitalize()+': '+strattag.capitalize(), fontsize=20)

ax.set_prop_cycle(color=[
    '#1f77b4', '#1f77b4', '#aec7e8', '#aec7e8', '#ff7f0e', '#ff7f0e', '#ffbb78', '#ffbb78', '#2ca02c', '#2ca02c', '#98df8a', '#98df8a',
    '#d62728', '#d62728', '#ff9896', '#ff9896', '#9467bd', '#9467bd', '#c5b0d5', '#c5b0d5'])#,
    # '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    # '#17becf', '#9edae5'])#,
    # marker=["s"]*20)


for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    # First test

    itest = tests[0]

    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
    track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'

    # Read variable
    varfil_main = Dataset(datdir+'strat.nc')
    strat = varfil_main.variables['strat'][:,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    varfil_main.close()
    nt = strat.shape[0]
    print(nt)
    sys.exit()

    var_ls = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)

    # Plot variable
    plt.plot(clon_shift, clat, linewidth=2, label=nustr[imemb])
    skip=24
    itim=np.arange(0,nt,skip)
    plt.plot(clon_shift[itim], clat[itim], "s")
    
    # Second test


plt.legend(loc="upper right")

# plt.show()
# plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'_'+memb_all[imemb]+'.png',dpi=200, facecolor='white', \
plt.savefig(figdir+'tser_'+storm+'_'+fig_extra+'.png',dpi=200, facecolor='white', \
            bbox_inches='tight', pad_inches=0.2)
plt.close()
