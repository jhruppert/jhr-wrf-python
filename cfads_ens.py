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
import sys
import cmocean
from thermo_functions import density_moist, theta_dry, theta_equiv, theta_virtual #, relh


# #### Variable selection

# Fill variable
iplot = 'vmf'
# options: vmf, thv, the

# Settings
# Calculate anomaly as deviation from xy-mean
do_prm_xy = 1
# Calculate anomaly as time-increment
do_prm_inc = 0

# Should be off for VMF
if iplot == 'vmf':
  do_prm_xy=0


# #### Test/storm selection

storm = 'haiyan'
#storm = 'maria'

nmem = 5 # number of ensemble members (1-5 have NCRF)
#nmem = 1

xmin=780


# #### Time selection

ntall=[1,3,6,12,24,36]
i_nt=np.shape(ntall)[0]

for knt in range(i_nt):
#for knt in range(3,4):
  
  nt = ntall[knt]
  hr_tag = str(np.char.zfill(str(nt), 2))
  
  
  # #### Directories
  
  figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'
  main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
  
  nums=np.arange(1,nmem+1,1); nums=nums.astype(str)
  nustr = np.char.zfill(nums, 2)
  memb_all=np.char.add('memb_',nustr)
  
  datdir2 = 'post/d02/'
  
  
  ##### Get dimensions
  
  datdir = main+storm+'/'+memb_all[0]+'/ctl/'+datdir2
  varfil_main = Dataset(datdir+'T.nc')
  tmpv = varfil_main.variables['T'][0:1,:,:,xmin:1400-1] # K
  pres = varfil_main.variables['pres'][:] # hPa
  varfil_main.close()
  bv_shape = np.shape(tmpv)
  nz = bv_shape[1]
  nx1 = bv_shape[2]
  nx2 = bv_shape[3]
  
  
  # Bin settings
  if iplot == 'thv':
    nbin=60
    fmin=-5; fmax=5
    step=(fmax-fmin)/nbin
    bins=np.arange(fmin,fmax,step)
  elif iplot == 'the':
    nbin=60
    fmin=-10; fmax=10
    step=(fmax-fmin)/nbin
    bins=np.arange(fmin,fmax,step)
  elif iplot == 'vmf':
    bins=np.logspace(-3,1.1,num=20)
    bins=np.concatenate((-1.*np.flip(bins),bins))
    nbin=np.shape(bins)[0]

  # Create axis of bin center-points
  bin_axis = (bins[np.arange(nbin-1)]+bins[np.arange(nbin-1)+1])/2


  # #### Read variables ##############################################
  
  ntest=2
  var_freq=np.zeros((ntest,nbin-1,nz))
  var_mn=np.zeros((ntest,nz))
  
  for ktest in range(ntest):
#  for ktest in range(1):
  
    if ktest == 0:
      itest='ctl'
      t0=36
    elif ktest == 1:
      itest='ncrf'
      t0=0
    t0=t0+1 # add one time step since NCRF(t=0) = CTL
    t1 = t0+nt
  
    print('Running itest: ',itest)
  
    # Create arrays to save ens members
    if do_prm_inc == 1:
      var_all = np.zeros((nmem,nt-1,nz,nx1,nx2))
    else:
      var_all = np.zeros((nmem,nt,nz,nx1,nx2))
  
    for imemb in range(nmem):
  
      print('Running imemb: ',memb_all[imemb])
  
      datdir = main+storm+'/'+memb_all[imemb]+'/'+itest+'/'+datdir2
  
  # Three-dimensional variables
  
  # Mixing ratio
      varfil_main = Dataset(datdir+'QVAPOR.nc')
      qv = varfil_main.variables['QVAPOR'][t0:t1,:,:,xmin:1400-1] # kg/kg
      varfil_main.close()
  
  # Temperature
      varfil_main = Dataset(datdir+'T.nc')
      tmpk = varfil_main.variables['T'][t0:t1,:,:,xmin:1400-1] # K
      varfil_main.close()
  
  
  ### Variable selection ##############################################
  
  # Virtual potential temp
      if iplot == 'thv':
        var = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
  
        # Figure settings
        fig_title=r"$\theta_v$"
        fig_tag='thv'
        units_var='K'
  
        # For mean var
        scale_mn=1.#e3
        units_mn=units_var
        xrange_mn=(-0.5,0.5)
        xrange_mn2=(-0.1,0.1)

  # Equiv potential temp
      elif iplot == 'the': 
        var = theta_equiv(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
  
        # Figure settings
        fig_title=r"$\theta_e$"
        fig_tag='the'
        units_var='K'
  
  # Vertical mass flux
      elif iplot == 'vmf':
        # Density
        rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3
        varfil_cvar = Dataset(datdir+'W.nc') # this opens the netcdf file
        var = varfil_cvar.variables['W'][t0:t1,:,:,xmin:1400-1] # m/s
        varfil_cvar.close()
        var *= rho
  
        # Figure settings
        fig_title='VMF'
        fig_tag='vmf'
        units_var='kg m$^{-2}$ s$^{-1}$'

        # For mean var
        scale_mn=1e3
        units_mn='$10^{-3}$ '+units_var
        xrange_mn=(-1,8)
        xrange_mn2=(-1,1)
  
  # Th_v' weighted by qv'
  #thv = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
  #thp = thv[range(1,nt),:,:,:] - thv[range(0,nt-1),:,:,:]
  #qvp = np.absolute(qv[range(1,nt),:,:,:] - qv[range(0,nt-1),:,:,:])
  #qvscale = np.mean(qvp,axis=(2,3))
  #qvp /= qvscale[:,:,np.newaxis,np.newaxis]
  #thp *= qvp
  
  # Calculate var' as time-increment: var[t,z,y,x] - mean_xy(var[t,z])
      if do_prm_xy == 1:
        v_mean = np.mean(var,axis=(2,3))
        var -= v_mean[:,:,np.newaxis,np.newaxis]
        fig_tag+='_xyp'
        fig_title+=' (xp)'
  
  # Calculate var' as time-increment: var[t] - var[t-1]
      if do_prm_inc == 1:
        tmpvar = var[range(1,nt),:,:,:] - var[range(0,nt-1),:,:,:]
        del var
        var = tmpvar
        fig_tag+='_tp'
        fig_title+=' (tp)'
  
  # Save ens member
      var_all[imemb,:,:,:,:] = var
  
  
  #### Calculate frequency ##############################################
  
    for iz in range(nz):
        for ibin in range(nbin-1):
            indices = ((var_all[:,:,iz,:,:] >= bins[ibin]) & (var_all[:,:,iz,:,:] < bins[ibin+1])).nonzero()
            var_freq[ktest,ibin,iz]=np.shape(indices)[1]
    
    ncell=nx1*nx2*nt*nmem
    var_freq[ktest,:,:] *= 100./ncell

  #### Calculate basic mean
    var_mn[ktest,:] = np.mean(var_all,axis=(0,1,3,4))  

  
  # ### Plotting routines ##############################################
  
  font = {'family' : 'sans-serif',
          'weight' : 'normal',
          'size'   : 16}
  
  matplotlib.rc('font', **font)
  
  
  # ### Plot CFADs for both tests ########################
  
  for ktest in range(ntest):
  
    if ktest == 0:
      itest='ctl'
    elif ktest == 1:
      itest='ncrf'

    pltvar = np.transpose(np.ma.masked_equal(var_freq[ktest,:,:],0))
    var_mn_plt = var_mn[ktest,:]*scale_mn

    fig, axd = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                            constrained_layout=True, figsize=(12, 8))

    ifig_title=fig_title+' ('+itest.upper()+')'
    fig.suptitle(ifig_title)

    for col in range(2):
        
        ax = plt.subplot(1,2,1+col)
        ax.set_yscale('log')
        ax.invert_yaxis()
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xlabel(units_var)
        ax.tick_params(axis='both',length=7)
        plt.yticks(ticks=pres)
        plt.ylim(np.max(pres), np.min(pres))


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
            elif iplot == 'thv':
                clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
                ticks=None
    
            var_plt=np.transpose(np.ma.masked_equal(var_freq,0))

            im = ax.contourf(bin_axis, pres, pltvar, clevs, norm=colors.LogNorm(),
                             cmap=cmocean.cm.ice_r, alpha=1.0, extend='max', zorder=2)
            
            plt.xlim(np.min(bin_axis), np.max(bin_axis))
            locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
            ax.xaxis.set_major_locator(locmin)
            # ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            
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

    plt.savefig(figdir+'cfad_'+fig_tag+'_ens5m_'+itest+'_'+hr_tag+'.png',dpi=200, facecolor='white', \
                bbox_inches='tight', pad_inches=0.2)

  
  # ### Plot difference CFAD ########################
  
  # var_diff = NCRF - CTL
  pltvar = np.transpose( var_freq[0,:,:] - var_freq[1,:,:] )
  var_mn_plt = (var_mn[0,:] - var_mn[1,:])*scale_mn
  
  fig, axd = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1]},
                          constrained_layout=True, figsize=(12, 8))
  
  ifig_title=fig_title+' (CTL - NCRF)'
  fig.suptitle(ifig_title)
  
  for col in range(2):
  
      ax = plt.subplot(1,2,1+col)
      ax.set_yscale('log')
      ax.invert_yaxis()
      ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
      ax.set_xlabel(units_var)
      ax.tick_params(axis='both',length=7)
      plt.yticks(ticks=pres)
      plt.ylim(np.max(pres), np.min(pres))
  
  
  ####### Fill contour ##############
  
      if col == 0:
  
          ax.set_title('CFAD')
          ax.set_ylabel('Pressure [hPa]')
          ax.set_xscale('symlog')
  
          clevsi=np.concatenate(([1e-2],np.arange(2,11,2)*1e-2,np.arange(2,11,2)*1e-1,np.arange(2,11,2)*1e-0))
          clevs = np.concatenate((-1*np.flip(clevsi),clevsi))
        
          im = ax.contourf(bin_axis, pres, pltvar, clevs, norm=colors.SymLogNorm(base=10,linthresh=clevsi[0],linscale=clevsi[0]), \
               cmap='RdBu_r', alpha=1, extend='max', zorder=2)
  
          plt.xlim(np.min(bin_axis), np.max(bin_axis))
          locmin = ticker.SymmetricalLogLocator(base=10.0,linthresh=2,subs=np.arange(2,11,2)*0.1)
          ax.xaxis.set_major_locator(locmin)
          # ax.xaxis.set_minor_formatter(ticker.NullFormatter())
  
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
  
  plt.savefig(figdir+'cfad_'+fig_tag+'_ens5m_diff_'+hr_tag+'.png',dpi=200, facecolor='white', \
              bbox_inches='tight', pad_inches=0.2)
  
