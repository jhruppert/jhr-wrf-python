#!/usr/bin/env python
# coding: utf-8

# ### Notebook to genereate histograms comparing tests across TC ens simulations
# Based on script cfads_ens_strat.py
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 4/23/22


from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker, cm
import subprocess
import sys
import cmocean
from thermo_functions import density_moist, theta_equiv, theta_virtual, relh
from mask_tc_track import mask_tc_track


# How many ensemble members
nmem = 10 # number of ensemble members (1-10 have NCRF)
# nmem = 5
# nmem = 1
enstag = str(nmem)
# Starting member to read
memb0=1

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
rmax = 8 # radius (deg) limit for masking around TC center


# #### Variable selection

# Fill variable
ivar_select = 'rain'
# options: pw, vmf, rain, lwacre

# #### Time selection

# ntall=[1,3,6,12,24,36]
# ntall=[1,3,6,12]
ntall=[1,6,12]
ntall=[1]

# #### Storm selection

# storm = 'haiyan'
# storm = 'maria'
storm_all=['haiyan','maria']
# storm_all=['haiyan']
# storm_all=['maria']
nstorm=np.size(storm_all)


print("Variable: ",ivar_select)


# istrat=2 # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
# for istrat in range(-1,4):
for istrat in range(2,3):
# for istrat in range(-1,3):

  print("Strat = ",istrat)

  # Index AKA Bin variable settings

  # PW
  if ivar_select == 'pw':
      fmin=35;fmax=80 # mm
      step=1
      bins=np.arange(fmin,fmax+step,step)
      xlabel='Column water vapor [mm]'
      log_x='linear'
  # Rainfall rate
  elif ivar_select == 'rain':
      bins=10.**(np.arange(1,8,0.3)-4)
      bins=10.**(np.arange(0,8,0.3)-4)
      xlabel='Rainfall rate [mm/hr]'
      log_x='log'
  # Vertical mass flux
  elif ivar_select == 'vmf':
      bins=10.**(np.arange(1,8,0.3)-3)
      # bins=np.flip(-1.*bins)
      xlabel='Vertical mass flux [kg/m/s]'
      log_x='log'
  # LW-ACRE
  elif ivar_select == 'lwacre':
      fmin=-50; fmax=200 # W/m2
      step=5
      bins=np.arange(fmin,fmax+step,step)
      xlabel='LW-ACRE [W/m**2]'
      log_x='linear'

  nbins = np.size(bins)
  # Create axis of bin center-points for plotting
  bin_axis = (bins[np.arange(nbins-1)]+bins[np.arange(nbins-1)+1])/2
  # Whereas "bins" holds the bin edges

  # Function for reading 2D variables
  def var_read_2d(datdir,varname,t0,t1):
      varfil_main = Dataset(datdir+varname+'.nc')
      var = varfil_main.variables[varname][t0:t1,:,:,:]
      varfil_main.close()
      return var


  for istorm in range(nstorm):

    storm=storm_all[istorm]
    print("Storm: ",storm)

    # Tests to compare
    if storm == 'haiyan':
      tests = ['ctl','ncrf36h']
    elif storm == 'maria':
      # tests = ['ctl','ncrf36h']
      tests = ['ctl','ncrf48h']
    # tests = ['crfon','ncrf']

    # Shift starting-read time step for CRFON comparison
    t0_test=0
    if tests[0] == 'crfon':
      t0_test=24
      memb0=5 # for CRFFON test

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
      print(strattag)

    # #### Time selection
    i_nt=np.shape(ntall)[0]

    for knt in range(i_nt):
    # for knt in range(0,1):
    #for knt in range(3,i_nt+1):

      nt = ntall[knt]
      hr_tag = str(np.char.zfill(str(nt), 2))
      print("Hour sample: ",hr_tag)

      # #### Directories

      figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'
      # main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
      main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

      nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
      nustr = np.char.zfill(nums, 2)
      memb_all=np.char.add('memb_',nustr)

      # datdir2 = 'post/d02/v2/'
      datdir2 = 'post/d02/'

      ##### Get dimensions

      datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
      varfil_main = Dataset(datdir+'T.nc')
      nz = varfil_main.dimensions['level'].size
      # lat = varfil_main.variables['XLAT'][:][0] # deg
      # lon = varfil_main.variables['XLONG'][:][0] # deg
      nx1 = varfil_main.dimensions['lat'].size
      nx2 = varfil_main.dimensions['lon'].size#-xmin-1
      pres = varfil_main.variables['pres'][:] # hPa
      varfil_main.close()

      process = subprocess.Popen(['ls '+main+storm+'/'+memb_all[0]+'/'+tests[0]+'/wrfout_d02_*'],shell=True,
          stdout=subprocess.PIPE,universal_newlines=True)
      output = process.stdout.readline()
      wrffil = output.strip() #[3]
      varfil_main = Dataset(wrffil)
      lat = varfil_main.variables['XLAT'][:][0] # deg
      lon = varfil_main.variables['XLONG'][:][0] # deg
      varfil_main.close()


      # #### Read variables ##############################################
      
      ntest=2
      var_freq=np.zeros((ntest,nbins-1,nz))
      
      for ktest in range(ntest):
      
        itest=tests[ktest]

        # if ktest == 0:
        #   nmem=1
        # else:
        #   nmem=5

        # This has been tested for corresponding time steps:
        #   t0=37,1 are the first divergent time steps in CTL,NCRF
        #   t0=25,1 are the first divergent time steps in NCRF,CRFON
        if itest == 'ctl':
          if tests[1] == 'ncrf36h':
            t0=36
          elif tests[1] == 'ncrf48h':
            t0=48
        elif itest == 'ncrf36h':
          t0=t0_test
        elif itest == 'ncrf48h':
          t0=t0_test
        elif itest == 'crfon':
          t0=0

        print('Running itest: ',itest)

        # Create arrays to save ens members
        var_all = np.ma.zeros((nmem,nt,nx1,nx2))
          # var_copy = np.ma.zeros((nmem,nt,nz,nx1,nx2))

        for imemb in range(nmem):
      
          print('Running imemb: ',memb_all[imemb])
      
          datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'+datdir2
          print(datdir)
      
          # Localize to TC track
          track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'

        # Two-dimensional variables

        # Stratiform index
          if istrat != -1:
            varfil_main = Dataset(datdir+'strat.nc')
            strat = varfil_main.variables['strat'][t0:t1,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
            varfil_main.close()      

        # Index AKA Bin variable ("ivar")

          # PW
          if ivar_select == 'pw':
              varname = ivar_select.upper()
              ivar = var_read_2d(datdir,varname,t0,t1)
          # Rainfall rate
          elif ivar_select == 'rain':
              ivar = rain # mm/hr
          # Vertical mass flux
          elif ivar_select == 'vmf':
              g=9.81 # gravity, m/s2
              wv_int = np.sum(w,axis=1) * dp/g # m/s * s**2/m * kg/m/s**2 = kg/s/m
              ivar = np.reshape(wv_int,(nt,1,nx1,nx2))
          # LW-ACRE
          elif ivar_select == 'lwacre':
              binfil = Dataset(datdir+'LWacre.nc') # this opens the netcdf file
              ivar = binfil.variables['LWUPB'][t0:t1,:,:,:] # W/m2
              binfil.close()

          # Mask out based on strat/conv
          if istrat != -1:
            var = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), var, copy=True)
            # vmf_copy = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), vmf_copy, copy=True)

          # Localize to TC track
          var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
          # vmf_copy = mask_tc_track(track_file, rmax, vmf_copy, lon, lat, t0, t1)

          # Save ens member
          var_all[imemb,:,:,:] = var
      
      #### Calculate frequency ##############################################




      # ### Plotting routines ##############################################
      
      font = {'family' : 'sans-serif',
              'weight' : 'normal',
              'size'   : 16}
      
      matplotlib.rc('font', **font)
      
      
      # ### Plot CFADs for both tests ########################
      
      for ktest in range(ntest):
      
        itest=tests[ktest]

        # pltvar = np.transpose(var_freq[ktest,:,:])
        pltvar = np.transpose(np.ma.masked_equal(var_freq[ktest,:,:],0))
        var_mn_plt = var_mn[ktest,:]*scale_mn

        fig, axd = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                                constrained_layout=True, figsize=(12, 8))

        ifig_title=fig_title+' ('+itest.upper()+')'
        if istrat != -1:
            ifig_title+=' ('+strattag+')'
        fig.suptitle(ifig_title)

        for col in range(2):
            
            ax = plt.subplot(1,2,1+col)

            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.tick_params(axis='both',length=7)
            ytick_loc=np.arange(900,0,-100)
            plt.yticks(ticks=ytick_loc)
            plt.ylim(np.max(pres), 100)#np.min(pres))

            ax.set_xlabel(units_var)


        ####### Fill contour ##############

            if col == 0:

                ax.set_title('CFAD')
                ax.set_ylabel('Pressure [hPa]')

                if iplot == 'vmf':
                    ax.set_xscale('symlog')
                    clevs=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)))
                    
                    locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
                    ax.xaxis.set_major_locator(locmin)
                    ticks=[1e-2,1e-1,1,1e1]
                else: #if iplot == 'thv' or iplot == 'the':
                    clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
                    ticks=None
        
                im = ax.contourf(bin_axis, pres, pltvar, clevs, norm=colors.LogNorm(),
                                cmap=cmocean.cm.ice_r, alpha=1.0, extend='max', zorder=2)
                
                plt.xlim(np.min(bin_axis), np.max(bin_axis))
                
                # ax2 = ax.twinx()
                # plt.ylim()
                # plt.plot(bin_axis,var_freq_int[ktest,:])

                cbar = plt.colorbar(im, ax=ax, shrink=0.75, ticks=ticks, format=ticker.LogFormatterMathtext())
                cbar.ax.set_ylabel('%')


        ####### Mean profile ##############

            elif col == 1:
        
                ax.set_title('Mean')
                ax.yaxis.set_major_formatter(ticker.NullFormatter())
                
                ax.plot(var_mn_plt, pres, "-k", linewidth=2)
                plt.xlim(xrange_mn)
                plt.axvline(x=0,color='k',linewidth=0.5)
                ax.set_xlabel(units_mn)

        plt.savefig(figdir+'cfad_'+fig_tag+fig_extra+'_ens'+enstag+'m_'+itest+'_'+hr_tag+'.png',dpi=200, facecolor='white', \
                    bbox_inches='tight', pad_inches=0.2)
        plt.close()



      
      # ### Plot difference CFAD ########################
      
      # var_diff = CTL - NCRF
      pltvar = np.transpose( var_freq[0,:,:] - var_freq[1,:,:] )
      var_mn_plt = (var_mn[0,:] - var_mn[1,:])*scale_mn
      
      fig, axd = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                              constrained_layout=True, figsize=(12, 8))
      
      ifig_title=fig_title+' ('+tests[0].upper()+' - '+tests[1].upper()+')'
      if istrat != -1:
          ifig_title+=' ('+strattag+')'
      fig.suptitle(ifig_title)
      
      for col in range(2):
      
          ax = plt.subplot(1,2,1+col)

          ax.set_yscale('log')
          ax.invert_yaxis()
          ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
          ax.tick_params(axis='both',length=7)
          ytick_loc=np.arange(900,0,-100)
          plt.yticks(ticks=ytick_loc)
          plt.ylim(np.max(pres), 100)#np.min(pres))

          ax.set_xlabel(units_var)
      
      
      ####### Fill contour ##############
      
          if col == 0:
      
              ax.set_title('CFAD')
              ax.set_ylabel('Pressure [hPa]')

              if iplot == 'vmf':
                  ax.set_xscale('symlog')
                  clevsi=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)*1e-0))

                  locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
                  ax.xaxis.set_major_locator(locmin)
              else: #if iplot == 'thv' or iplot == 'the':
                  if iplot == 'qrad':
                    clevsi=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)*1e0,np.arange(2,11,2)*1e1))
                  else:
                    clevsi=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)*1e0))

              clevs = np.concatenate((-1*np.flip(clevsi),clevsi))
              ticks=[1e-2,1e-1,1,1e1]

              im = ax.contourf(bin_axis, pres, pltvar, clevs, norm=colors.SymLogNorm(base=10,linthresh=clevsi[0],linscale=clevsi[0]),
                              cmap='RdBu_r', alpha=1.0, extend='max', zorder=2)

              plt.xlim(np.min(bin_axis), np.max(bin_axis))

              #if iplot == 'thv': 
              plt.axvline(x=0,color='k',linewidth=1.)

              cbar = plt.colorbar(im, ax=ax, shrink=0.75, ticks=ticker.SymmetricalLogLocator(base=10.0, linthresh=.5),
                                  format=ticker.LogFormatterMathtext())
              cbar.ax.set_ylabel('%')

        ####### Mean profile ##############

          elif col == 1:

              ax.set_title('Mean')
              ax.yaxis.set_major_formatter(ticker.NullFormatter())

              ax.plot(var_mn_plt, pres, "-k", linewidth=2)
              plt.xlim(xrange_mn2)
              plt.axvline(x=0,color='k',linewidth=0.5)
              ax.set_xlabel(units_mn)

      difftag='diff'
      if tests[0] == 'crfon': difftag+='v2'
      plt.savefig(figdir+'cfad_'+fig_tag+fig_extra+'_ens'+enstag+'m_'+difftag+'_'+hr_tag+'.png',dpi=200, facecolor='white', \
                  bbox_inches='tight', pad_inches=0.2)
      plt.close()
