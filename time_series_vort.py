#!/usr/bin/env python
# coding: utf-8

# ### Notebook to plot time series of vorticity/max wind
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/8/23

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
# storms=['maria']
# storms=['haiyan']
# storm = 'haiyan'
# storm = 'maria'

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

    # How many members
    nmem = 10 # number of ensemble members
    # nmem = 2
    # Starting member to read

    # Strat/Conv index subset
#     istrat=1 # Convective
# #    istrat=2 # Stratiform
#     istrat=4 # Convective/Stratiform fraction
#     istrat_all=[1,2,4]
#     nstrat=np.size(istrat_all)


    # TC tracking
    ptrack=600 # tracking pressure level
    var_track = 'rvor' # variable
    # rmax = 6 # radius (deg) limit for masking around TC center

    # #### Directories

    figdir = "/home/jamesrup/figures/tc/ens/vorticity/"
    main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

    memb0=1
    nums=np.arange(memb0,nmem+memb0,1)
    nums=nums.astype(str)
    nustr = np.char.zfill(nums, 2)
    memb_all=np.char.add('memb_',nustr)

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

    # Level selection
    ikread = np.where(pres == ptrack)[0][0]

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

    nt = np.zeros(ntest, dtype=np.int32)
    for itest in range(ntest):
        ##### Get dimensions
        datdir = main+storm+'/'+memb_all[0]+'/'+tests[itest]+'/'
        varfil_main = Dataset(datdir+'post/d02/T.nc')
        i_nt = varfil_main.dimensions['time'].size
        varfil_main.close()
        nt[itest]=i_nt


    mlvort_t0 = np.zeros((nmem,nt[0]))
    vmax_t0 = np.zeros((nmem,nt[0]))
    satfrac_t0 = np.zeros((nmem,nt[0]))

    mlvort_t1 = np.zeros((nmem,nt[1]))
    vmax_t1 = np.zeros((nmem,nt[1]))
    satfrac_t1 = np.zeros((nmem,nt[1]))

    for itest in range(ntest):
    # for itest in range(1):

        print('Running test: ',tests[itest])

        for imemb in range(nmem):

            print('Running imemb: ',memb_all[imemb])

            # First test

            datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/'

            # track_file = datdir+'track_'+var_track+'_'+ptrack+'hPa.nc'
            # Localize to TC track
            # NOTE: Using copied tracking from CTL for NCRF tests
            trackfil_ex=''
            if 'ncrf' in tests[itest]:
                trackfil_ex='_ctlcopy'
            track_file = datdir+'track_'+var_track+trackfil_ex+'_'+str(round(pres[ikread]))+'hPa.nc'

            # Read variables

            datdir = main+storm+'/'+memb_all[imemb]+'/'+tests[itest]+'/post/d02/'

            # Winds
            varfil_main = Dataset(datdir+'U10.nc')
            u10 = varfil_main.variables['U10'][:,:,:,:] # m/s
            varfil_main.close()
            varfil_main = Dataset(datdir+'V10.nc')
            v10 = varfil_main.variables['V10'][:,:,:,:] # m/s
            varfil_main.close()
            wsp = np.sqrt(u10**2 + v10**2)

            # Midlevel vorticity
            fil = Dataset(datdir+'AVOR.nc') # this opens the netcdf file
            vort = fil.variables['AVOR'][:,ikread,:,:] # 10**-5 /s
            fil.close()
            shape_sav=vort.shape
            vort = np.reshape(vort, (shape_sav[0],1,shape_sav[1],shape_sav[2]))

            # Saturation fraction
            varfil_main = Dataset(datdir+'satfrac.nc')
            pw = varfil_main.variables['pw'][:,:,:,:] # mm
            pws = varfil_main.variables['pw_sat'][:,:,:,:] # mm
            varfil_main.close()
            satfraci = pw/pws

            # Mask out around TC center
            t0_test=0
            t1_test=nt[itest]
            rmax = 2 # radius (deg) limit for masking around TC center
            vort = mask_tc_track(track_file, rmax, vort, lon, lat, t0_test, t1_test)
            wsp = mask_tc_track(track_file, rmax, wsp, lon, lat, t0_test, t1_test)
            satfraci = mask_tc_track(track_file, rmax, satfraci, lon, lat, t0_test, t1_test)

            # Average / take max

            vortmax = np.mean(vort, axis=(2,3))
            vortmax = np.reshape(vortmax,nt[itest])
            vortmax = np.ma.filled(vortmax, fill_value=np.nan)

            wsp = np.ma.filled(wsp, fill_value=np.nan)
            wspmax = np.nanmax(wsp, axis=(2,3))
            wspmax = np.reshape(wspmax,nt[itest])

            satfrac = np.mean(satfraci, axis=(2,3))
            satfrac = np.reshape(satfrac,nt[itest])
            satfrac = np.ma.filled(satfrac, fill_value=np.nan)

            if itest == 0:
                mlvort_t0[imemb,:] = vortmax
                vmax_t0[imemb,:] = wspmax
                satfrac_t0[imemb,:] = satfrac
            elif itest == 1:
                mlvort_t1[imemb,:] = vortmax
                vmax_t1[imemb,:] = wspmax
                satfrac_t1[imemb,:] = satfrac


    # ### Plotting routines ##############################################

    font = {'family' : 'sans-serif',
            'weight' : 'normal',
            'size'   : 14}

    matplotlib.rc('font', **font)

    nvar=3
    for ivar in range(nvar):

        if ivar == 0:
            var0 = np.copy(mlvort_t0)
            var1 = np.copy(mlvort_t1)
            title_tag = 'Midlevel Vorticity'
            figtag = 'mlvort'
            ylabel = '10$^{-5}$ /s'
        elif ivar == 1:
            var0 = np.copy(vmax_t0)
            var1 = np.copy(vmax_t1)
            title_tag = 'Max 10m Wind'
            figtag = 'vmax'
            ylabel = 'm/s'
        elif ivar == 2:
            var0 = np.copy(satfrac_t0)
            var1 = np.copy(satfrac_t1)
            title_tag = 'Saturation fraction'
            figtag = 'satfrac'
            ylabel = '%'

        var0 = pd.DataFrame(var0)
        var0 = var0.rolling(window=3, center=True, closed='both', axis=1).mean()
        var0 = np.copy(var0)

        var1 = pd.DataFrame(var1)
        var1 = var1.rolling(window=3, center=True, closed='both', axis=1).mean()
        var1 = np.copy(var1)

        # create figure
        fig = plt.figure(figsize=(9,5))
        ax = fig.add_subplot(111)

        ax.set_title(title_tag+' ('+storm.capitalize()+')')#, fontsize=20)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Time [hours]')

        t_range=[0,96]
        plt.xlim(t_range)

        color_t0 = 'red'
        color_t1 = 'blue'

        # Test 0

        mean_t0 = np.nanmean(var0, axis=0)
        std_t0 = np.nanstd(var0, axis=0)

        tshift = get_tshift(tests[0])
        xdim = range(0+tshift, nt[0]+tshift)

        # for imemb in range(nmem):
        #     plt.plot(xdim, var0[imemb,:], linewidth=2, label=tests[0].upper(), color=color_t0, linestyle='solid')
        plt.plot(xdim, mean_t0, linewidth=2, label=tests[0].upper(), color=color_t0, linestyle='solid')
        plt.fill_between(xdim, mean_t0 + std_t0, mean_t0 - std_t0, alpha=0.2, color=color_t0)

        # Test 1

        mean_t1 = np.nanmean(var1, axis=0)
        std_t1 = np.nanstd(var1, axis=0)

        tshift = get_tshift(tests[1])
        xdim = range(0+tshift, nt[1]+tshift)

        # for imemb in range(nmem):
        #     plt.plot(xdim, var1[imemb,:], linewidth=2, label=tests[0].upper(), color=color_t1, linestyle='solid')
        plt.plot(xdim, mean_t1, linewidth=2, label=tests[0].upper(), color=color_t1, linestyle='solid')
        plt.fill_between(xdim, mean_t1 + std_t1, mean_t1 - std_t1, alpha=0.2, color=color_t1)

        plt.grid()

        # plt.legend(loc="upper right")

        plt.savefig(figdir+'tser_'+storm+'_'+figtag+'.png',dpi=200, facecolor='white', \
                    bbox_inches='tight', pad_inches=0.2)
        # plt.show()
        plt.close()