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
from matplotlib import ticker, cm
import sys
#import cmocean
from thermo_functions import density_moist, theta_dry, theta_equiv, theta_virtual #, relh
# from stratiform_index import stratiform_index


# #### Variable selection

# Fill variable
iplot = 'vmf'
# options: vmf, thv, the

# Settings
# Calculate anomaly as deviation from xy-mean
do_prm_xy = 0
# Calculate anomaly as time-increment
do_prm_inc = 0

# Should be off for VMF
if iplot == 'vmf':
  do_prm_xy=0


# #### Test/storm selection

storm = 'haiyan'
#storm = 'maria'

nmem = 2 #5 # number of ensemble members (1-5 have NCRF)


# #### Time selection

nt=6

xmin=780


# #### Directories

#figdir = "/Users/jruppert/code/tc_figs/"
figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'

# main = "/Users/jamesruppert/code/tc_output/"
# main = "/Users/jruppert/code/tc_output/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"

nums=np.arange(1,nmem+1,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)
#nmem=np.shape(memb_all)[0]
#imemb=memb_all[memb-1]
#memb = get_ipython().getoutput('ls $main/$istorm')
#imemb=memb[0]
# print(main+istorm+'/'+imemb)

#datdir = main+storm+'/'+imemb+'/'+itest+'/'
datdir2 = 'post/d02/'


# Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/ctl/'+datdir2
varfil_main = Dataset(datdir+'T.nc')
tmpv = varfil_main.variables['T'][0:1,:,:,xmin:1400-1] # K
pres = varfil_main.variables['pres'][:] # hPa
varfil_main.close()
bv_shape = np.shape(tmpv)
nz = bv_shape[1]
nx1 = bv_shape[2]
nx2 = bv_shape[3]


# #### Read variables ##############################################

ntest=2
nbin=60
var_freq=np.zeros((ntest,nbin-1,nz))

for ktest in range(ntest):

  if ktest == 0:
    itest='ctl'
    t0=36
  elif ktest == 1:
    itest='ncrf'
    t0=0
  t1 = t0+nt

  print('Running itest: ',itest)

  # Create arrays to save ens members
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

      # CFAD settings
      fmin=-5; fmax=5

      # Figure settings
      fig_title="i$\theta_v$"
      fig_tag='thv'
      units_var='K'

# Equiv potential temp
    elif iplot == 'the': 
      var = theta_equiv(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

      # CFAD settings
      fmin=-10; fmax=10

      # Figure settings
      fig_title="i$\theta_e$"
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

      # CFAD settings
      fmin=-3; fmax=5
      if do_prm_inc == 1:
        fmin=-5; fmax=5

      # Figure settings
      fig_title="VMF'"
      fig_tag='vmf'
      units_var='kg m$^{-2}$ s$^{-1}$'

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

# Calculate var' as time-increment: var[t] - var[t-1]
    if do_prm_inc == 1:
      tmpvar = var[range(1,nt),:,:,:] - var[range(0,nt-1),:,:,:]
      del var
      var = tmpvar
      fig_tag+='_tp'

# Save ens member
    var_all[imemb,:,:,:,:] = var


#### Calculate frequency ##############################################

  step=(fmax-fmin)/nbin
  bins=np.arange(fmin,fmax,step)
 
  for iz in range(nz):
      for ibin in range(nbin-1):
          indices = ((var_all[:,:,iz,:,:] >= bins[ibin]) & (var_all[:,:,iz,:,:] < bins[ibin+1])).nonzero()
          var_freq[ktest,ibin,iz]=np.shape(indices)[1]
  
  ncell=nx1*nx2*nt*nmem
  # print(ncell)
  var_freq[ktest,:,:] *= 100./ncell



# ### Plotting routines ##############################################

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


# ### Plot CFADs for both tests

for ktest in range(ntest):

  if ktest == 0:
    itest='ctl'
  elif ktest == 1:
    itest='ncrf'

  pltvar = np.transpose(var_freq[ktest,:,:])

  fig = plt.figure(figsize=(14,8))
  ax = fig.add_subplot(111)
  
  ax.set_title(fig_title)
  ax.set_ylabel('Pressure [hPa]')
  
  ax.invert_yaxis()
  ax.set_yscale('log')
  
  # fill contour
  clevs=[0.005,0.01,0.05,0.1,0.5,1,5,10,50]
  
  im = ax.contourf(bins[0:nbin-1], pres, pltvar, clevs, locator=ticker.LogLocator(), \
                   cmap='GnBu', alpha=1, extend='max', zorder=2)
  
  ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
  
  plt.xlim(np.min(bins), np.max(bins))
  
  plt.axvline(x=0,color='k',linewidth=0.5)
  ax.set_xlabel(units_var)
  
  cbar = plt.colorbar(im, ax=ax, shrink=0.75)
  cbar.ax.set_ylabel('%')
  
  plt.savefig(figdir+'cfad_'+fig_tag+'_ens5m_'+itest+'.png',dpi=200, facecolor='white', \
              bbox_inches='tight', pad_inches=0.2)



# ### Plot difference CFAD

# var_diff = NCRF - CTL
pltvar = np.transpose( var_freq[1,:,:] - var_freq[0,:,:] )

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111)

ax.set_title(fig_title)
ax.set_ylabel('Pressure [hPa]')

ax.invert_yaxis()
ax.set_yscale('log')

# fill contour
clevs=[0.005,0.01,0.05,0.1,0.5,1,5,10,50]
clevs = np.concatenate((-1*np.flip(clevs),clevs))

im = ax.contourf(bins[0:nbin-1], pres, pltvar, clevs, locator=ticker.LogLocator(), \
                 cmap='GnBu', alpha=1, extend='max', zorder=2)

ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())

plt.xlim(np.min(bins), np.max(bins))

plt.axvline(x=0,color='k',linewidth=0.5)
ax.set_xlabel(units_var)

cbar = plt.colorbar(im, ax=ax, shrink=0.75)
cbar.ax.set_ylabel('%')

plt.savefig(figdir+'cfad_'+fig_tag+'_ens5m_diff.png',dpi=200, facecolor='white', \
            bbox_inches='tight', pad_inches=0.2)


