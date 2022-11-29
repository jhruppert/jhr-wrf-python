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
storm = 'maria'

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
elif storm == 'maria':
    tests = ['ctl','ncrf48h']

# How many members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 2
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

figdir = "/home/jamesrup/figures/tc/ens/"
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

def get_tshift(itest):
    if itest == 'ctl':
        tshift=0
    elif itest == 'ncrf36h':
        tshift=36
    elif itest == 'ncrf48h':
        tshift=48
    return tshift


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
ax.set_ylabel('Fraction')
ax.set_xlabel('Time [hours]')

# ax.set_prop_cycle(color=[
#     '#1f77b4', '#1f77b4', 
#     '#aec7e8', '#aec7e8', 
#     '#ff7f0e', '#ff7f0e', 
#     '#ffbb78', '#ffbb78', 
#     '#2ca02c', '#2ca02c', 
#     '#98df8a', '#98df8a',
#     '#d62728', '#d62728', 
#     '#ff9896', '#ff9896', 
#     '#9467bd', '#9467bd', 
#     '#c5b0d5', '#c5b0d5',
#     '#000000', '#000000'])#,
#     # '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
#     # '#17becf', '#9edae5'])#,
#     # marker=["s"]*20)

color_t1 = 'red'
color_t2 = 'blue'

for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    # First test

    itest = tests[0]
    tshift1 = get_tshift(itest)

    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'
    track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'

    # Read variable
    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
    varfil_main = Dataset(datdir+'strat.nc')
    strat = varfil_main.variables['strat'][:,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    varfil_main.close()
    nt = strat.shape[0]

    t0_test1=0
    t1_test1=nt

    # Mask out around TC center
    strat = mask_tc_track(track_file, rmax, strat, lon, lat, t0_test1, t1_test1)
    count_total = np.ma.MaskedArray.count(strat, axis=(1,2,3))

    # Mask out based on strat/conv
    strat = np.ma.masked_where((strat != istrat), strat, copy=True)
    count_strat = np.ma.MaskedArray.count(strat, axis=(1,2,3))
    frac_strat = count_strat / count_total

    # Plot variable
    plt.plot(range(t0_test1+tshift1,t1_test1+tshift1), frac_strat, linewidth=1, 
        label=nustr[imemb], color=color_t1, linestyle='solid')

    if imemb == 0:
        frac_strat_all_t1 = np.reshape(frac_strat, (nt,1))
    else:
        frac_strat_all_t1 = np.append(frac_strat_all_t1, np.reshape(frac_strat, (nt,1)), axis=1)


    # Second test

    itest = tests[1]
    tshift2 = get_tshift(itest)

    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'
    track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'

    # Read variable
    datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
    varfil_main = Dataset(datdir+'strat.nc')
    strat = varfil_main.variables['strat'][:,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    varfil_main.close()
    nt = strat.shape[0]

    t0_test2=0
    t1_test2=nt

    # Mask out around TC center
    strat = mask_tc_track(track_file, rmax, strat, lon, lat, t0_test2, t1_test2)
    count_total = np.ma.MaskedArray.count(strat, axis=(1,2,3))

    # Mask out based on strat/conv
    strat = np.ma.masked_where((strat != istrat), strat, copy=True)
    count_strat = np.ma.MaskedArray.count(strat, axis=(1,2,3))
    frac_strat = count_strat / count_total

    # Plot variable
    plt.plot(range(t0_test2+tshift2,t1_test2+tshift2), frac_strat, linewidth=1,
        label=nustr[imemb], color=color_t2, linestyle='--')

    if imemb == 0:
        frac_strat_all_t2 = np.reshape(frac_strat, (nt,1))
    else:
        frac_strat_all_t2 = np.append(frac_strat_all_t2, np.reshape(frac_strat, (nt,1)), axis=1)

# Plot means

frac_mean_t1 = np.mean(frac_strat_all_t1, axis=1)
plt.plot(range(t0_test1 + tshift1, t1_test1 + tshift1), frac_mean_t1, 
    linewidth=2, label=nustr[imemb], color=color_t1, linestyle='solid')

frac_mean_t2 = np.mean(frac_strat_all_t2, axis=1)
plt.plot(range(t0_test2 + tshift2, t1_test2 + tshift2), frac_mean_t2, 
    linewidth=2, label=nustr[imemb], color=color_t2, linestyle='--')

# plt.legend(loc="upper right")

# plt.show()
# plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'_'+memb_all[imemb]+'.png',dpi=200, facecolor='white', \
plt.savefig(figdir+'tser_'+storm+fig_extra+'.png',dpi=200, facecolor='white', \
            bbox_inches='tight', pad_inches=0.2)
plt.close()
