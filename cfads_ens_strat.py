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
from mask_tc_track import mask_tc_track
from cfads_functions import cfads_var_settings, cfads_var_calc


# How many ensemble members
nmem = 10 # number of ensemble members (1-10 have NCRF)
# nmem = 3
enstag = str(nmem)

# #### Variable selection

# Fill variable
ivar_all = ['thv','vmf','lh','rh','qrad']
ivar_all = ['thv','vmf','lh','rh']
# ivar_all = ['lh','rh']
ivar_all = ['wpthp','wpthep','vmf','thv','the']
ivar_all = ['wpthp','wpthep']
ivar_all = ['lh','thv','the']
ivar_all = ['lh','vmf']
ivar_all = ['lq']
ivar_all = ['lh','vmf']
ivar_all = ['vmf']
nvar=np.size(ivar_all)

# #### Time selection

# ntall=[1,3,6,12,24,36]
# ntall=[1,3,6,12]
# ntall=[1,6,12]
ntall=[1,2,3,6]
ntall=[1,6]
ntall=[1]
# ntall=[1,3,6]

# #### Classification selection

# 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
# kclass=[0,1,2,3]
kclass=[2]

# #### Storm selection

# storm = 'haiyan'
# storm = 'maria'
storm_all=['haiyan','maria']
storm_all=['haiyan']
# storm_all=['maria']
nstorm=np.size(storm_all)

# TC tracking
ptrack='600' # tracking pressure level
var_track = 'rvor' # variable
# rmax = 8 # radius (deg) limit for masking around TC center
rmax = 3


for ivar in range(nvar):
# for ivar in range(1):

  iplot = ivar_all[ivar]
  print("Variable: ",iplot)

  # Calculate anomaly as deviation from xy-mean
  do_prm_xy = 0
  # Calculate anomaly as time-increment
  do_prm_inc = 0
  # if (iplot == 'thv') or (iplot == 'qrad'):
  if (iplot == 'thv') or (iplot == 'the') or (iplot == 'lq'):
      do_prm_xy = 1
  # Should be off for VMF
  if (iplot == 'vmf') or ('wpth' in iplot):
      do_prm_xy=0


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

      memb0=1
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

      # Get variable settings
      bins, fig_title, fig_tag, units_var, scale_mn, \
        units_mn, xrange_mn, xrange_mn2 = cfads_var_settings(iplot)
      nbin=np.shape(bins)[0]

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

      # Create arrays to save ens members
      if do_prm_inc == 1:
        var_all = np.ma.zeros((ntest,nmem,nt+1,nz,nx1,nx2)) # time dim will be reduced to nt in the subtraction
      else:
        var_all = np.ma.zeros((ntest,nmem,nt,nz,nx1,nx2))
        # var_copy = np.ma.zeros((ntest,nmem,nt,nz,nx1,nx2))
        strat_all = np.ma.zeros((ntest,nmem,nt,1,nx1,nx2))

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

        for imemb in range(nmem):
      
          print('Running imemb: ',memb_all[imemb])
      
          datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'+datdir2
          print(datdir)
      
          # Localize to TC track
          # track_file = datdir+'../../../track_'+var_track+'_'+ptrack+'hPa.nc'
          # track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'
          # Localize to TC track
          # NOTE: Using copied tracking from CTL for NCRF tests
          # track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'
          trackfil_ex=''
          if 'ncrf' in itest:
              trackfil_ex='_ctlcopy'
          track_file = datdir+'../../track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        # Two-dimensional variables

        # Stratiform index
          # if istrat != -1:
          varfil_main = Dataset(datdir+'strat.nc')
          strat = varfil_main.variables['strat'][t0:t1,:,:,:] # 0-non-raining, 1-conv, 2-strat, 3-other/anvil
          varfil_main.close()

        # Three-dimensional variables
          var = cfads_var_calc(iplot, datdir, pres, t0, t1)

        ### Process variable ##############################################

          # Calculate var' as anomaly from x-y-average, using large-scale (large-radius) var avg
          if do_prm_xy == 1:
            radius_ls=6
            var_ls = mask_tc_track(track_file, radius_ls, var, lon, lat, t0, t1)
            # var_ls = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
            var_ls_avg = np.ma.mean(var_ls,axis=(0,2,3))
            var -= var_ls_avg[np.newaxis,:,np.newaxis,np.newaxis]

          # Calculate w'th'
          if ('wpth' in iplot):
            print("NOT RUNNING THIS VAR. NEED TO FIX READ FUNCTION")
            print("PUNTING TO FUTURE SELF WHO IS WISER THAN I")
            sys.exit()
            # rmax_prime = rmax
            # var_ls1 = mask_tc_track(track_file, rmax_prime, thp, lon, lat, t0, t1)
            # var_ls1_avg = np.ma.mean(var_ls1,axis=(0,2,3))
            # thp -= var_ls1_avg[np.newaxis,:,np.newaxis,np.newaxis]
            # var_ls2 = mask_tc_track(track_file, rmax_prime, www, lon, lat, t0, t1)
            # var_ls2_avg = np.ma.mean(var_ls2,axis=(0,2,3))
            # www -= var_ls2_avg[np.newaxis,:,np.newaxis,np.newaxis]
            # var = thp*www

          # # Mask out based on strat/conv
          # if (istrat != -1) & (istrat != 3):
          #   var = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), var, copy=True)
          #   # vmf_copy = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), vmf_copy, copy=True)
          # elif istrat == 3:
          #   var = np.ma.masked_where(((np.repeat(strat,nz,axis=1) == 0) & (np.repeat(strat,nz,axis=1) == 3)),
          #       var, copy=True)

          # Localize to TC track
          var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
          # vmf_copy = mask_tc_track(track_file, rmax, vmf_copy, lon, lat, t0, t1)

          # Save ens member
          var_all[ktest,imemb,:,:,:,:] = var
          strat_all[ktest,imemb,:,:,:,:] = strat
          # var_copy[imemb,:,:,:,:] = vmf_copy

      # Calculate var' as time-increment: var[t] - var[t-1]
      if do_prm_inc == 1:
        var_all = var_all[:,:,1:,:,:,:] - var_all[:,:,:-1,:,:,:]

##### CLASSIFICATION LOOP ###############

      # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
      nclass=np.shape(kclass)[0]
      for kstrat in range(nclass):
        
        istrat=kclass[kstrat]

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
            # strattag='Anv'
            strattag='CS'
          fig_extra='_'+strattag.lower()
          print("Strat tag: ",strattag)

        # Mask out based on strat/conv
        if ((istrat != -1) & (istrat != 3)):
          var_ma = np.ma.masked_where((np.repeat(strat_all,nz,axis=3) != istrat), var_all, copy=True)
          # vmf_copy = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), vmf_copy, copy=True)
        elif istrat == 3:
          var_ma = np.ma.masked_where(((np.repeat(strat_all,nz,axis=3) == 0) & (np.repeat(strat_all,nz,axis=3) == 3)),
              var_all, copy=True)

        #### Calculate basic mean
        var_mn=np.ma.mean(var_ma, axis=(1,2,4,5))

#### Calculate frequency ##############################################

        for ktest in range(ntest):
          for iz in range(nz):
            for ibin in range(nbin-1):
              indices = ((var_ma[ktest,:,:,iz,:,:] >= bins[ibin]) & (var_ma[ktest,:,:,iz,:,:] < bins[ibin+1])).nonzero()
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
        
        
# ### Plot CFADs for both tests ################################
        
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

                  if (iplot == 'vmf') or ('wpth' in iplot):
                      ax.set_xscale('symlog')
                      clevs=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)))
                      
                      locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
                      # if iplot == 'vmf':
                      #   locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
                      # elif iplot == 'wpthp':
                      #   locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
                      ax.xaxis.set_major_locator(locmin)
                      ticks=[1e-2,1e-1,1,1e1]
                  else: #if iplot == 'thv' or iplot == 'the':
                      clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
                      # clevs=[1,2,5,10,25,50,100,200,500]
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

                if (iplot == 'vmf') or ('wpth' in iplot):
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
                # ticks=[1e-2,1e-1,1,1e1]

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
