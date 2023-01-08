#!/usr/bin/env python
# coding: utf-8

# #### Notebook to create time series of precipitation efficiency (PE) from TC output
# 
# Assumes key integrated variables have been processed and written out by PE_write.py to single-level netcdf files.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/5/23


from netCDF4 import Dataset
import numpy as np
import matplotlib
# matplotlib.use('pdf')
import matplotlib.pyplot as plt
import subprocess
import sys
from mask_tc_track import mask_tc_track
import pandas as pd


# #### Main settings


# NOTE: Using copied tracking from CTL for NCRF tests

# #### Variable selection

storm = 'haiyan'
storm = 'maria'

# How many members
nmem = 10 # number of ensemble members
# nmem = 4

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 6 # radius (deg) limit for masking around TC center

# Strat/Conv index subset
# istrat_all=[0,1,2] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
nrain=5 # np.size(istrat_all)
    # krain = 
    # 0 = conv+strat points
    # 1 = conv points
    # 2 = strat points
    # 3 = rainfall rate threshold
    # 4 = all points (inside of TC mask)

# #### Directories

figdir = "/home/jamesrup/figures/tc/ens/precip_eff/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"


def get_tshift(itest):
    if itest == 'ctl':
        tshift=0
    elif itest == 'ncrf36h':
        tshift=36
    elif itest == 'ncrf48h':
        tshift=48
    return tshift


# Tests to read and compare
ntest=2
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
elif storm == 'maria':
#        tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']

# Ens member strings
memb0=1 # Starting member to read
nums=np.arange(memb0,nmem+memb0,1)
nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

# Get Lat/Lon
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
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


# Create arrays

nt = np.zeros(ntest, dtype=np.int32)

for itest in range(ntest):
    ##### Get dimensions
    datdir = main+storm+'/'+memb_all[0]+'/'+tests[itest]+'/'
    varfil_main = Dataset(datdir+'post/d02/T.nc')
    i_nt = varfil_main.dimensions['time'].size
    varfil_main.close()

    nt[itest]=i_nt

mf_ratio_t0 = np.zeros((nmem,nrain,nt[0]))
pe_mf_t0 = np.zeros((nmem,nrain,nt[0]))
pe_mp_t0 = np.zeros((nmem,nrain,nt[0]))

mf_ratio_t1 = np.zeros((nmem,nrain,nt[1]))
pe_mf_t1 = np.zeros((nmem,nrain,nt[1]))
pe_mp_t1 = np.zeros((nmem,nrain,nt[1]))


def plot_rainhist(x):
    n, bins, patches = plt.hist(x, 500, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('mm/hr')
    plt.ylabel('Occurrence')
    plt.title('Rainfall Rate Distribution')
    plt.xlim(0.1, 80)
    # plt.ylim(0, 0.03)
    # plt.grid(True)
    # plt.show()


# #### Main loop

for itest in range(ntest):

    print('Running test: ',tests[itest])

    for imemb in range(nmem):

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/'
        # track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
        # Localize to TC track
        # NOTE: Using copied tracking from CTL for NCRF tests
        trackfil_ex=''
        if 'ncrf' in tests[itest]:
            trackfil_ex='_ctlcopy'
        track_file = datdir+'track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        # Read variables

        # Strat
        datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/post/d02/'
        varfil_main = Dataset(datdir+'strat.nc')
        strat = varfil_main.variables['strat'][:,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
        varfil_main.close()

        # Rain
        datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/post/d02/'
        varfil_main = Dataset(datdir+'rainrate.nc')
        rain = varfil_main.variables['rainrate'][:,:,:,:]/24 # mm/d --> mm/hr
        varfil_main.close()

        # rain = np.ma.masked_where(rain <= 3, rain)
        # print(rain.min())
        # plot_rainhist(rain.flatten())
        # sys.exit()

        # PE variables
        datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/post/d02/'
        varfil_main = Dataset(datdir+'precip_eff_vars.nc')
        vmfu = varfil_main.variables['vmfu'][:,:,:,:] # kg/m/s
        vmfd = varfil_main.variables['vmfd'][:,:,:,:] # kg/m/s
        condh = varfil_main.variables['condh'][:,:,:,:] # mm/d
        condh /= 24 # mm/d --> mm/hr
        varfil_main.close()

        t0=0
        t1=nt[itest]

        # Mask out around TC center
        rain = mask_tc_track(track_file, rmax, rain, lon, lat, t0, t1)
        rain = np.ma.filled(rain, fill_value=np.nan)
        strat = mask_tc_track(track_file, rmax, strat, lon, lat, t0, t1)
        strat = np.ma.filled(strat, fill_value=np.nan)

        # PE terms
        mf_ratio = -1 * vmfd / vmfu
        pe_mf = 1 - mf_ratio
        pe_mp = rain / condh
        mf_ratio = np.ma.masked_invalid(mf_ratio)
        pe_mf = np.ma.masked_invalid(pe_mf)
        pe_mp = np.ma.masked_invalid(pe_mp)

        # Average across raining points
        for it in range(nt[itest]):
            for krain in range(nrain):
                # krain=3

                # conv+strat points
                if krain == 0:
                    ind_rain = ((strat[it,0,:,:] == 1) | (strat[it,0,:,:] == 2)).nonzero()
                # conv points
                elif krain == 1:
                    ind_rain = (strat[it,0,:,:] == 1).nonzero()
                # strat points
                elif krain == 2:
                    ind_rain = (strat[it,0,:,:] == 2).nonzero()
                # rainfall rate threshold
                elif krain == 3:
                    rain_thresh = 3. # mm/hr
                    ind_rain = (rain[it,0,:,:] >= rain_thresh).nonzero()
                # Simple mean across (unmasked) domain
                # elif krain == 4:

                # if krain < 4:
                #     condh_avg = np.nanmean(condh[it,0,ind_rain[0],ind_rain[1]])
                #     rain_avg = np.nanmean(rain[it,0,ind_rain[0],ind_rain[1]])
                #     vmfu_avg = np.nanmean(vmfu[it,0,ind_rain[0],ind_rain[1]])
                #     vmfd_avg = np.nanmean(vmfd[it,0,ind_rain[0],ind_rain[1]])
                # else:
                #     condh_avg = np.nanmean(condh[it,0,:,:])
                #     rain_avg = np.nanmean(rain[it,0,:,:])
                #     vmfu_avg = np.nanmean(vmfu[it,0,:,:])
                #     vmfd_avg = np.nanmean(vmfd[it,0,:,:])

                # mf_ratio = -1 * vmfd_avg / vmfu_avg
                # pe_mf = 1 - mf_ratio
                # pe_mp = rain_avg / condh_avg

                # if itest == 0:
                #     mf_ratio_t0[imemb,krain,it] = mf_ratio
                #     pe_mf_t0[imemb,krain,it] = pe_mf
                #     pe_mp_t0[imemb,krain,it] = pe_mp
                # elif itest == 1:
                #     mf_ratio_t1[imemb,krain,it] = mf_ratio
                #     pe_mf_t1[imemb,krain,it] = pe_mf
                #     pe_mp_t1[imemb,krain,it] = pe_mp

                if krain < 4:
                    mf_ratio_avg = np.mean(mf_ratio[it,0,ind_rain[0],ind_rain[1]])
                    pe_mf_avg = np.mean(pe_mf[it,0,ind_rain[0],ind_rain[1]])
                    pe_mp_avg = np.mean(pe_mp[it,0,ind_rain[0],ind_rain[1]])
                else:
                    mf_ratio_avg = np.mean(mf_ratio[it,0,:,:])
                    pe_mf_avg = np.mean(pe_mf[it,0,:,:])
                    pe_mp_avg = np.mean(pe_mp[it,0,:,:])

                if itest == 0:
                    mf_ratio_t0[imemb,krain,it] = mf_ratio_avg
                    pe_mf_t0[imemb,krain,it] = pe_mf_avg
                    pe_mp_t0[imemb,krain,it] = pe_mp_avg
                elif itest == 1:
                    mf_ratio_t1[imemb,krain,it] = mf_ratio_avg
                    pe_mf_t1[imemb,krain,it] = pe_mf_avg
                    pe_mp_t1[imemb,krain,it] = pe_mp_avg


# ---
# ### Plotting routines


font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)


for krain in range(nrain):
# for krain in range(1):

    # conv+strat points
    if krain == 0:
        fig_extra='CS'
        raintag='Conv|Strat'
    # conv points
    elif krain == 1:
        fig_extra='conv'
        raintag='Convective'
    # strat points
    elif krain == 2:
        fig_extra='strat'
        raintag='Stratiform'
    # rainfall rate threshold
    elif krain == 3:
        fig_extra='rainthresh'
        raintag='Rainfall threshold'
    # rainfall rate threshold
    elif krain == 4:
        fig_extra='all'
        raintag='All points'

    mf0 = mf_ratio_t0[:,krain,:]
    pe_mf0 = pe_mf_t0[:,krain,:]
    pe_mp0 = pe_mp_t0[:,krain,:]

    mf1 = mf_ratio_t1[:,krain,:]
    pe_mf1 = pe_mf_t1[:,krain,:]
    pe_mp1 = pe_mp_t1[:,krain,:]

    nvar=3
    for ivar in range(nvar):
    # for ivar in range(0,1):

        if ivar == 0:
            var0 = mf0
            var1 = mf1
            title_tag = 'MF Ratio (dn/up)'
            figtag = 'mffrac'
        elif ivar == 1:
            var0 = pe_mf0
            var1 = pe_mf1
            title_tag = 'PE (MF)'
            figtag = 'pemf'
        elif ivar == 2:
            var0 = pe_mp0
            var1 = pe_mp1
            title_tag = 'PE (MP)'
            figtag = 'pemp'

        # var0 = pd.DataFrame(var0)
        # # var0 = 
        # var0.rolling(window=3, center=True, closed='both', axis=0).mean()
        # var1 = pd.DataFrame(var1)
        # # var1 = 
        # var1.rolling(window=3, center=True, closed='both', axis=0).mean()

        # create figure
        fig = plt.figure(figsize=(9,5))
        ax = fig.add_subplot(111)

        ax.set_title(title_tag+' ('+storm.capitalize()+'; '+raintag+')')#, fontsize=20)
        ax.set_ylabel('Fraction')
        ax.set_xlabel('Time [hours]')

        t_range=[30,80]
        # plt.xlim(t_range)

        color_t0 = 'red'
        color_t1 = 'blue'

        # Test 0

        mean_t0 = np.nanmean(var0, axis=0)
        std_t0 = np.nanstd(var0, axis=0)

        tshift = get_tshift(tests[0])
        xdim = range(0+tshift, nt[0]+tshift)

        plt.plot(xdim, mean_t0, linewidth=2, label=tests[0].upper(), color=color_t0, linestyle='solid')
        plt.fill_between(xdim, mean_t0 + std_t0, mean_t0 - std_t0, alpha=0.2, color=color_t0)

        # Test 1

        mean_t1 = np.nanmean(var1, axis=0)
        std_t1 = np.nanstd(var1, axis=0)

        tshift = get_tshift(tests[1])
        xdim = range(0+tshift, nt[1]+tshift)

        plt.plot(xdim, mean_t1, linewidth=2, label=tests[0].upper(), color=color_t1, linestyle='solid')
        plt.fill_between(xdim, mean_t1 + std_t1, mean_t1 - std_t1, alpha=0.2, color=color_t1)
        plt.grid()

        # plt.legend(loc="upper right")

        plt.savefig(figdir+'tser_'+storm+'_'+figtag+'_'+fig_extra+'.png',dpi=200, facecolor='white',                     bbox_inches='tight', pad_inches=0.2)
        # plt.show()
        plt.close()

