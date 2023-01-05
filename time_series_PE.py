#!/usr/bin/env python
# coding: utf-8

# ### Notebook to plot time series of precipitation efficiency quantified in different ways.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/4/23

# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
from matplotlib import ticker, cm
import subprocess
from mask_tc_track import mask_tc_track
import sys
import pandas as pd


# #### Variable selection

storms=['haiyan','maria']
storms=['maria']
storms=['haiyan']
# storm = 'haiyan'
# storm = 'maria'

# How many members
nmem = 10 # number of ensemble members
# nmem = 1
# Starting member to read
memb0=1

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 6 # radius (deg) limit for masking around TC center

# Strat/Conv index subset
# istrat_all=[0,1,2] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
nrain=4 # np.size(istrat_all)
    # irain = 
    # 0 = conv+strat points
    # 1 = conv points
    # 2 = strat points
    # 3 = rainfall rate threshold


def get_strattag(istrat):
    # Strat/Conv index subset
    if istrat == -1:
        strattag='all'
    elif istrat == 0:
        strattag='nonrain'
    elif istrat == 1:
        strattag='conv'
    elif istrat == 2:
        strattag='strat'
    elif istrat == 3:
        strattag='anv'

def get_tshift(itest):
    if itest == 'ctl':
        tshift=0
    elif itest == 'ncrf36h':
        tshift=36
    elif itest == 'ncrf48h':
        tshift=48
    return tshift


nstorm = np.size(storms)
for istorm in range(nstorm):

    storm = storms[istorm]

    # Tests to read and compare
    ntest=2
    if storm == 'haiyan':
        tests = ['ctl','ncrf36h']
    elif storm == 'maria':
#        tests = ['ctl','ncrf36h']
        tests = ['ctl','ncrf48h']

    # #### Directories

    figdir = "/home/jamesrup/figures/tc/ens/"
    main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

    nums=np.arange(memb0,nmem+memb0,1)
    nums=nums.astype(str)
    nustr = np.char.zfill(nums, 2)
    memb_all=np.char.add('memb_',nustr)

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


    mf_ratio = np.zeros((ntest,nmem,nrain,nt))
    pw_mf = np.zeros((ntest,nmem,nrain,nt))
    pe_mp = np.zeros((ntest,nmem,nrain,nt))

    for imemb in range(nmem):

        print('Running imemb: ',memb_all[imemb])

        for itest in range(ntest):

            itest = tests[itest]
            tshift1 = get_tshift(itest)

            datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'
            # track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
            # Localize to TC track
            # NOTE: Using copied tracking from CTL for NCRF tests
            trackfil_ex=''
            if 'ncrf' in itest:
                trackfil_ex='_ctlcopy'
            track_file = datdir+'track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

            # Read variables

            # Strat
            datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
            varfil_main = Dataset(datdir+'strat.nc')
            strat = varfil_main.variables['strat'][:,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
            varfil_main.close()
            nt = strat.shape[0]

            # Rain
            datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
            varfil_main = Dataset(datdir+'rainrate.nc')
            rain = varfil_main.variables['rainrate'][:,:,:,:] # mm/d
            varfil_main.close()

            # PE variables
            datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/'
            varfil_main = Dataset(datdir+'precip_eff_vars.nc')
            vmfu = varfil_main.variables['vmfu'][:,:,:,:] # kg/m/s
            vmfd = varfil_main.variables['vmfd'][:,:,:,:] # kg/m/s
            condh = varfil_main.variables['condh'][:,:,:,:] # mm/d
            varfil_main.close()

            t0_test1=0
            t1_test1=nt

            # Mask out around TC center
            rain = mask_tc_track(track_file, rmax, rain, lon, lat, t0_test1, t1_test1)
            strat = mask_tc_track(track_file, rmax, strat, lon, lat, t0_test1, t1_test1)

            # Average across raining points
            for it in range(nt):
                for krain in range(nrain):

                    # conv+strat points
                    if krain == 0:
                        irain = (np.logical_or((strat[it,0,:,:] == 1) , (strat[it,0,:,:] == 2))).nonzero()
                    # conv points
                    elif krain == 1:
                        irain = (strat[it,0,:,:] == 1).nonzero()
                    # strat points
                    elif krain == 2:
                        irain = (strat[it,0,:,:] == 2).nonzero()
                    # rainfall rate threshold
                    elif krain == 3:
                        rain_thresh = 1. # mm/d
                        irain = (rain[it,0,:,:] >= rain_thresh).nonzero()

                    vmfu_avg = np.mean(vmfu[it,0,irain[0],irain[1]])
                    vmfd_avg = np.mean(vmfd[it,0,irain[0],irain[1]])
                    condh_avg = np.mean(condh[it,0,irain[0],irain[1]])
                    rain_avg = np.mean(rain[it,0,irain[0],irain[1]])

                    mf_ratio[itest,imemb,krain,it] = vmfd_avg / vmfu_avg
                    pw_mf[itest,imemb,krain,it] = 1 - mf_ratio
                    pe_mp[itest,imemb,krain,it] = rain_avg / condh_avg

        sys.exit()


    # ### Plotting routines ##############################################

    # Create figure with all tracks

    for iplot in range(3):

        if iplot == 0:
            fig_extra='strat'
            strattag='Stratiform'
            pvar1 = frac_strat_all_t1
            pvar2 = frac_strat_all_t2
        elif iplot == 1:
            fig_extra='conv'
            strattag='Convective'
            pvar1 = frac_conv_all_t1
            pvar2 = frac_conv_all_t2
        elif iplot == 2:
            fig_extra='ratio'
            strattag='Convective/Stratiform'
            pvar1 = frac_conv_all_t1 / frac_strat_all_t1
            pvar2 = frac_conv_all_t2 / frac_strat_all_t2

        pvar_pd1 = pd.DataFrame(pvar1)
        pvar1_smooth = pvar_pd1.rolling(window=3, center=True, closed='both', axis=0).mean()
        pvar_pd2 = pd.DataFrame(pvar2)
        pvar2_smooth = pvar_pd2.rolling(window=3, center=True, closed='both', axis=0).mean()

        font = {'family' : 'sans-serif',
                'weight' : 'normal',
                'size'   : 16}

        matplotlib.rc('font', **font)


        # ### Combined plot: STRAT ##############################################

        # create figure
        fig = plt.figure(figsize=(9,5))
        ax = fig.add_subplot(111)

        t_range=[30,80]

        ax.set_title(storm.capitalize()+': '+strattag, fontsize=20)
        ax.set_ylabel('Fraction')
        ax.set_xlabel('Time [hours]')
        plt.xlim(t_range)

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

        # Plot means

        frac_mean_t1 = np.nanmean(pvar1_smooth, axis=1)
        frac_std_t1 = np.nanstd(pvar1_smooth, axis=1)

        plt.plot(range(t0_test1 + tshift1, t1_test1 + tshift1), frac_mean_t1, 
            linewidth=2, label=tests[0].upper(), color=color_t1, linestyle='solid')
        plt.fill_between(range(t0_test1 + tshift1, t1_test1 + tshift1), frac_mean_t1 + frac_std_t1,
            frac_mean_t1 - frac_std_t1, alpha=0.2, color=color_t1)


        frac_mean_t2 = np.nanmean(pvar2_smooth, axis=1)
        frac_std_t2 = np.nanstd(pvar2_smooth, axis=1)

        plt.plot(range(t0_test2 + tshift2, t1_test2 + tshift2), frac_mean_t2, 
            linewidth=2, label=tests[1].upper(), color=color_t2, linestyle='--')
        plt.fill_between(range(t0_test2 + tshift2, t1_test2 + tshift2), frac_mean_t2 + frac_std_t2,
            frac_mean_t2 - frac_std_t2, alpha=0.2, color=color_t2)

        # plt.legend(loc="upper right")

        # plt.show()
        # plt.savefig(figdir+storm+'_track_'+var_track+'_'+ptrack+'_'+memb_all[imemb]+'.png',dpi=200, facecolor='white', \
        plt.savefig(figdir+'tser_'+storm+'_'+fig_extra+'_'+tests[1]+'.png',dpi=200, facecolor='white', \
                    bbox_inches='tight', pad_inches=0.2)
        plt.close()
