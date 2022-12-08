#!/usr/bin/env python
# coding: utf-8

# ### Notebook to genereate CFADS from TC output
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
iplot = 'thv'#'the'#'vmf'#'rh'#'qrad'#
iplot = 'vmf'
iplot = 'lh'

ivar_all = ['thv','vmf','lh','rh','qrad']
ivar_all = ['thv','vmf','lh','rh']
ivar_all = ['thv']
nvar=np.size(ivar_all)

# #### Time selection

# ntall=[1,3,6,12,24,36]
# ntall=[1,3,6,12]
ntall=[1,6,12]
ntall=[1]

# #### Storm selection

# storm = 'haiyan'
# storm = 'maria'
storm_all=['haiyan','maria']
storm_all=['haiyan']
storm_all=['maria']
nstorm=np.size(storm_all)


for ivar in range(nvar):
# for ivar in range(1):

  iplot = ivar_all[ivar]
  print("Variable: ",iplot)

  # Calculate anomaly as deviation from xy-mean
  do_prm_xy = 0
  # Calculate anomaly as time-increment
  do_prm_inc = 0
  # if (iplot == 'thv') or (iplot == 'qrad'):
  if (iplot == 'thv'):
      do_prm_xy = 1
  # Should be off for VMF
  if iplot == 'vmf':
      do_prm_xy=0


  # istrat=2 # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
  # for istrat in range(-1,4):
  for istrat in range(2,3):
  # for istrat in range(-1,3):

    print("Strat = ",istrat)
    # continue


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

        # Variable settings

        if iplot == 'thv':

            # Bin settings
            nbin=60
            fmax=6 #5 #; fmin=-5
            #step=(fmax-fmin)/nbin
            step=fmax*2/nbin
            bins=np.arange(0,fmax,step)+step
            bins=np.concatenate((-1.*np.flip(bins),bins))
            nbin=np.shape(bins)[0]
        
            # Figure settings
            fig_title=r"$\theta_v$"
            fig_tag='thv'
            units_var='K'
        
            # For mean var
            scale_mn=1.
            units_mn=units_var
            xrange_mn=(-1,1)
            xrange_mn2=(-1,1)

        elif iplot == 'the':

            # Bin settings
            nbin=60
            fmax=15 #; fmin=-10
            #step=(fmax-fmin)/nbin
            step=fmax*2/nbin
            bins=np.arange(0,fmax,step)+step
            bins=np.concatenate((-1.*np.flip(bins),bins))
            nbin=np.shape(bins)[0]
        
            # Figure settings
            fig_title=r"$\theta_e$"
            fig_tag='the'
            units_var='K'
        
            # For mean var
            scale_mn=1.#e3
            units_mn=units_var
            xrange_mn=(-0.5,0.5)
            xrange_mn2=(-0.1,0.1)

        elif iplot == 'vmf':

            # Bin settings
            bins=np.logspace(-3,1.1,num=20)
            # bins=np.logspace(-3.5,0.7,num=20)
            bins=np.concatenate((-1.*np.flip(bins),bins))
            nbin=np.shape(bins)[0]

            # Figure settings
            fig_title='VMF'
            fig_tag='vmf'
            units_var='kg m$^{-2}$ s$^{-1}$'
            # fig_title='$w$'
            # fig_tag='w'
            # units_var='m s$^{-1}$'

            # For mean var
            scale_mn=1e2
            units_mn='$10^{-2}$ '+units_var
            xrange_mn=(-20,20)
            xrange_mn2=(-5,5)

        elif iplot == 'rh':

            # Bin settings
            # bins=np.logspace(-3,1.1,num=20)
            # bins=np.concatenate((-1.*np.flip(bins),bins))
            # nbin=np.shape(bins)[0]
            nbin=45
            fmax=125; fmin=-10
            step=(fmax-fmin)/nbin
            bins=np.arange(fmin,fmax,step)
            nbin=np.shape(bins)[0]

            # Figure settings
            fig_title='RH'
            fig_tag='rh'
            units_var='%'

            # For mean var
            scale_mn=1
            units_mn=units_var
            xrange_mn=(-1,105)
            xrange_mn2=(-2,2)
          
        elif iplot == 'qrad':

            # Bin settings
            nbin=60
            fmax=20 #; fmin=-10
            #step=(fmax-fmin)/nbin
            step=fmax*2/nbin
            bins=np.arange(0,fmax,step)+step
            bins=np.concatenate((-1.*np.flip(bins),bins))
            nbin=np.shape(bins)[0]

            # Figure settings
            fig_title='$Q_R$'
            fig_tag='qrad'
            units_var='K d$^{-1}$'

            # For mean var
            scale_mn=1.
            units_mn=units_var
            xrange_mn=(-8,3)
            xrange_mn2=(-10,7)

        elif iplot == 'lh':

            # Bin settings
            nbin=60
            fmax=20 #; fmin=-10
            #step=(fmax-fmin)/nbin
            step=fmax*2/nbin
            bins=np.arange(0,fmax,step)+step
            bins=np.concatenate((-1.*np.flip(bins),bins))
            nbin=np.shape(bins)[0]

            # Figure settings
            fig_title='$Q_L$'
            fig_tag='lheat'
            units_var='K hr$^{-1}$'

            # For mean var
            scale_mn=1.
            units_mn=units_var
            xrange_mn=(-1,1)
            xrange_mn2=(-0.5,0.5)

        # Create axis of bin center-points
        bin_axis = (bins[np.arange(nbin-1)]+bins[np.arange(nbin-1)+1])/2

        if do_prm_xy == 1:
            fig_tag+='_xyp'
            fig_title+="$'$"# (xp)'
        if do_prm_inc == 1:
            fig_tag+='_tp'
            fig_title+=' (tp)'


        # #### Read variables ##############################################
        
        ntest=2
        var_freq=np.ma.zeros((ntest,nbin-1,nz))
        var_freq_int=np.ma.zeros((ntest,nbin-1))
        var_mn=np.ma.zeros((ntest,nz))
        
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

          if do_prm_inc == 0:
            t0+=1 # add one time step since NCRF(t=0) = CTL

          t1 = t0+nt
          if do_prm_inc == 1:
            t1+=1

          print('Running itest: ',itest)

          # Create arrays to save ens members
          if do_prm_inc == 1:
            var_all = np.ma.zeros((nmem,nt+1,nz,nx1,nx2)) # time dim will be reduced to nt in the subtraction
          else:
            var_all = np.ma.zeros((nmem,nt,nz,nx1,nx2))
            # var_copy = np.ma.zeros((nmem,nt,nz,nx1,nx2))

          for imemb in range(nmem):
        
            print('Running imemb: ',memb_all[imemb])
        
            datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'+datdir2
            print(datdir)
        
            # Localize to TC track
            # track_file = datdir+'../../../track_'+var_track+'_'+ptrack+'hPa.nc'
            track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'

          # Two-dimensional variables

          # Stratiform index
            if istrat != -1:
              varfil_main = Dataset(datdir+'strat.nc')
              strat = varfil_main.variables['strat'][t0:t1,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
              varfil_main.close()
              # Testing strat/conv classification mods
              # varfil_main = Dataset(main+storm+'/'+memb_all[imemb]+'/'+itest+'/post/d02/v2/strat_origit1.nc') # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
              # strat1 = varfil_main.variables['strat'][:,:,:,:]
              # varfil_main.close()
              # varfil_main = Dataset(main+storm+'/'+memb_all[imemb]+'/'+itest+'/'+'post/d02/v2/strat.nc') # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
              # strat2 = varfil_main.variables['strat'][:,:,:,:]
              # varfil_main.close()
              # strat=strat2

          # Three-dimensional variables
        
          # Mixing ratio
            varfil_main = Dataset(datdir+'QVAPOR.nc')
            qv = varfil_main.variables['QVAPOR'][t0:t1,:,:,:] # kg/kg
            varfil_main.close()

          # Temperature
            varfil_main = Dataset(datdir+'T.nc')
            tmpk = varfil_main.variables['T'][t0:t1,:,:,:] # K
            varfil_main.close()
            
            
          ### Variable selection ##############################################
        
            # Virtual potential temp
            if iplot == 'thv':
              var = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
            # Equiv potential temp
            elif iplot == 'the': 
              var = theta_equiv(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
            # Vertical mass flux
            elif iplot == 'vmf':
              # Density
              rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3
              varfil = Dataset(datdir+'W.nc') # this opens the netcdf file
              var = varfil.variables['W'][t0:t1,:,:,:] # m/s
              varfil.close()
              # vmf_copy=np.copy(var)
              var *= rho
            # Humidity
            elif iplot == 'rh':
              var = relh(qv,pres[np.newaxis,:,np.newaxis,np.newaxis]*1e2,tmpk,ice=1) # %
            # Radiation
            elif iplot == 'qrad':
              varfil = Dataset(datdir+'RTHRATLW.nc') # this opens the netcdf file
              var = varfil.variables['RTHRATLW'][t0:t1,:,:,:]*3600*24 # K/s --> K/d
              varfil.close()
              varfil = Dataset(datdir+'RTHRATSW.nc') # this opens the netcdf file
              var += varfil.variables['RTHRATSW'][t0:t1,:,:,:]*3600*24 # K/s --> K/d
              varfil.close()
            # Latent heat
            elif iplot == 'lh':
              varfil = Dataset(datdir+'H_DIABATIC.nc') # this opens the netcdf file
              var = varfil.variables['H_DIABATIC'][t0:t1,:,:,:]*3600 # K/s --> K/hr
              varfil.close()

            ### Process variable ##############################################

            # Calculate var' as anomaly from x-y-average, using large-scale (large-radius) var avg
            if do_prm_xy == 1:
              # radius_ls=12
              # var_ls = mask_tc_track(track_file, radius_ls, var, lon, lat, t0, t1)
              var_ls = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
              var_ls_avg = np.ma.mean(var_ls,axis=(0,2,3))
              var -= var_ls_avg[np.newaxis,:,np.newaxis,np.newaxis]

            # Mask out based on strat/conv
            if istrat != -1:
              var = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), var, copy=True)
              # vmf_copy = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), vmf_copy, copy=True)

            # Localize to TC track
            var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
            # vmf_copy = mask_tc_track(track_file, rmax, vmf_copy, lon, lat, t0, t1)

            # Save ens member
            var_all[imemb,:,:,:,:] = var
            # var_copy[imemb,:,:,:,:] = vmf_copy

        #### Calculate basic mean
          var_mn[ktest,:]=np.ma.mean(var_all,axis=(0,1,3,4))

        # Calculate var' as time-increment: var[t] - var[t-1]
          if do_prm_inc == 1:
            var_all = var_all[:,1:,:,:,:] - var_all[:,:-1,:,:,:]
        
        
        #### Calculate frequency ##############################################

          for iz in range(nz):
            for ibin in range(nbin-1):
              indices = ((var_all[:,:,iz,:,:] >= bins[ibin]) & (var_all[:,:,iz,:,:] < bins[ibin+1])).nonzero()
              var_freq[ktest,ibin,iz]=np.shape(indices)[1]
            var_freq[ktest,:,iz] /= np.ma.sum(var_freq[ktest,:,iz])
          #ncell=nx1*nx2*nt*nmem
          var_freq[ktest,:,:] *= 100. #/ncell

          # Integrate binvar
          # dp = (pres[1]-pres[0])*1e2 # Pa
          # # p_int = [900,500] # layer to integrate over
          # p_int = [800,600] # layer to integrate over
          # ik0 = np.where(pres == p_int[0])[0][0]; ik1 = np.where(pres == p_int[1])[0][0]
          # var_copy.shape
          # var_int = np.sum(var_copy[:,:,ik0:ik1,:,:],axis=2) * dp / 9.81

          # Vertically integrated variable
          # for ibin in range(nbin-1):
          #   indices = ((var_int >= bins[ibin]) & (var_int < bins[ibin+1])).nonzero()
          #   var_freq_int[ktest,ibin]=np.shape(indices)[1]
          #   # print(np.shape(indices)[1])
          # var_freq_int[ktest,:] /= np.ma.sum(var_freq_int[ktest,:])
          # var_freq_int[ktest,:] *= 100.


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
