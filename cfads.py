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
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import sys
#import cmocean
from thermo_functions import density_moist, theta_dry, theta_equiv, theta_virtual, relh
# from stratiform_index import stratiform_index


# #### Variable selection

# Fill variable
fillvar_select = 'tprm'
# options: lwcrf, tprm, dbz, rh


# #### Test/storm selection

storm = 'haiyan'
#storm = 'maria'

memb=1 # 1-20
#memb=3

itest='ctl'
itest='ncrf' #'ctl'


# #### Time selection

nt=12
t0 = 0 #36 #48
t1 = t0+nt

xmin=780


# #### Directories

figdir = "/Users/jruppert/code/tc_figs/"
figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'

# main = "/Users/jamesruppert/code/tc_output/"
# main = "/Users/jruppert/code/tc_output/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"

#istorm=storm[0]
nums=np.arange(1,21,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)
imemb=memb_all[memb-1]
#memb = get_ipython().getoutput('ls $main/$istorm')
#imemb=memb[0]
# print(main+istorm+'/'+imemb)

datdir = main+storm+'/'+imemb+'/'+itest+'/'
datdir = datdir+'post/d02/'
print(datdir)


# #### Read variables

# Two-dimensional variables

# Conv/strat separation: varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain
#varfil_strat = Dataset(datdir+'strat.nc') # this opens the netcdf file
#strat = varfil_strat.variables['strat'][t0:t1,:,:,xmin:1400-1]
#varfil_strat.close()

## PW
#binfil = Dataset(datdir+'PW.nc') # this opens the netcdf file
#pw = binfil.variables['PW'][t0:t1,:,:,xmin:1400-1]
#binfil.close()
#
## LW-ACRE
#binfil = Dataset(datdir+'LWacre.nc') # this opens the netcdf file
#lwacre = binfil.variables['LWUPB'][t0:t1,:,:,xmin:1400-1] # W/m2
#binfil.close()
#
## Rainfall
#binfil = Dataset(datdir+'rainrate.nc') # this opens the netcdf file
#rain = binfil.variables['rainrate'][t0:t1,:,:,xmin:1400-1] # mm/hr
#binfil.close()



# Three-dimensional variables

# Vertical coordinate
filtmp = Dataset(datdir+'RTHRATLW.nc')
pres = filtmp.variables['pres'][:] # hPa
print("Vertical shape: ",np.shape(pres))
filtmp.close()

# Mixing ratio
varfil_main = Dataset(datdir+'QVAPOR.nc')
qv = varfil_main.variables['QVAPOR'][t0:t1,:,:,xmin:1400-1] # kg/kg
varfil_main.close()
units_var1 = 'kg/kg'

# Horizontal temperature anomaly
varfil_main = Dataset(datdir+'T.nc')
tmpk = varfil_main.variables['T'][t0:t1,:,:,xmin:1400-1] # K
varfil_main.close()
# thp = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
units_thp = 'K'
# Subtract time-dependent domain average
# t_mean = np.mean(np.mean(thp,axis=3),axis=2)
# thp -= t_mean[:,:,np.newaxis,np.newaxis]

# Equiv pot temp
# th_e = theta_equiv(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
# Subtract time-dependent domain average
# t_mean = np.mean(np.mean(th_e,axis=3),axis=2)
# th_e -= t_mean[:,:,np.newaxis,np.newaxis]

# Density
rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3

bv_shape = np.shape(tmpk)
print("Binvar shape: ",bv_shape)
nt = bv_shape[0]
nz = bv_shape[1]
nx1 = bv_shape[2]
nx2 = bv_shape[3]

# # Calculate th_v' as th_v[t] - th_v[t-1]
thv = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
thp = thv[range(1,nt),:,:,:] - thv[range(0,nt-1),:,:,:]

# # Calculate th_e' as th_e[t] - th_e[t-1]
# the = theta_equiv(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
# thp = the[range(1,nt),:,:,:] - the[range(0,nt-1),:,:,:]

# Th_v' weighted by qv'
#thv = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
#thp = thv[range(1,nt),:,:,:] - thv[range(0,nt-1),:,:,:]
#qvp = np.absolute(qv[range(1,nt),:,:,:] - qv[range(0,nt-1),:,:,:])
#qvscale = np.mean(qvp,axis=(2,3))
#qvp /= qvscale[:,:,np.newaxis,np.newaxis]
#thp *= qvp

# Vertical motion
varfil_cvar = Dataset(datdir+'W.nc') # this opens the netcdf file
w = varfil_cvar.variables['W'][t0:t1,:,:,xmin:1400-1] # m/s
varfil_cvar.close()
w *= rho
units_var2='kg m$^{-2}$ s$^{-1}$'
# lcmin = -20; lcmax=20; lcint=2

# # Calculate w' as w[t] - w[t-1]
#wp = w[range(1,nt),:,:,:] - w[range(0,nt-1),:,:,:]



# #### Calculate frequency

nbin=60

fmin=-3; fmax=5
#fmin=-5; fmax=5
step=(fmax-fmin)/nbin
wbins=np.arange(fmin,fmax,step)
w_freq=np.zeros((nbin-1,nz))

fmin=-5; fmax=5
#fmin=-10; fmax=10
step=(fmax-fmin)/nbin
thvbins=np.arange(fmin,fmax,step)
thv_freq=np.zeros((nbin-1,nz))

for iz in range(nz):
    for ibin in range(nbin-1):
        indices = ((w[:,iz,:,:] >= wbins[ibin]) & (w[:,iz,:,:] < wbins[ibin+1])).nonzero()
#        indices = ((wp[:,iz,:,:] >= wbins[ibin]) & (wp[:,iz,:,:] < wbins[ibin+1])).nonzero()
        w_freq[ibin,iz]=np.shape(indices)[1]
        # indices = ((thp[:,iz,:,:] >= thvbins[ibin]) & (thp[:,iz,:,:] < thvbins[ibin+1])).nonzero()
        # thv_freq[ibin,iz]=np.shape(indices)[1]

ncell=nx1*nx2*nt
# print(ncell)
w_freq *= 100./ncell
thv_freq *= 100./ncell


# ---
# ### Plotting routines

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


# ### Plot profiles

# In[13]:


# create figure
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111)

ax.set_title('VMF')
ax.set_ylabel('Pressure [hPa]')

ax.invert_yaxis()
ax.set_yscale('log')

# fill contour
clevs=[0.005,0.01,0.05,0.1,0.5,1,5,10,50]

im = ax.contourf(wbins[0:nbin-1], pres, np.transpose(w_freq), clevs, locator=ticker.LogLocator(), \
                 cmap='GnBu', alpha=1, extend='max', zorder=2)

# X=wbins[0:nbin-1]
# Y=pres
# Z=np.transpose(w_freq)
# ax.pcolor(X, Y, Z, norm=colors.LogNorm(vmin=clevs[0], vmax=np.max(clevs)),
#                    cmap='PuBu', shading='auto')

ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())

plt.xlim(np.min(wbins), np.max(wbins))

# plt.axvline(x=0,color='k',linewidth=0.5)
ax.set_xlabel(units_var2)

cbar = plt.colorbar(im, ax=ax, shrink=0.75)
cbar.ax.set_ylabel('%')

plt.savefig(figdir+'cfad_vmf_'+imemb+'_'+itest+'.png',dpi=200, facecolor='white', \
            bbox_inches='tight', pad_inches=0.2)


# In[84]:

#
## create figure
#fig = plt.figure(figsize=(14,8))
#ax = fig.add_subplot(111)
#
#ax.set_title("VMF'")
#ax.set_ylabel('Pressure [hPa]')
#
#ax.invert_yaxis()
#ax.set_yscale('log')
#
## fill contour
#clevs=[0.005,0.01,0.05,0.1,0.5,1,5,10,50]
#
#im = ax.contourf(wbins[0:nbin-1], pres, np.transpose(w_freq), clevs, locator=ticker.LogLocator(), \
#                 cmap='GnBu', alpha=1, extend='max', zorder=2)
#
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#
#plt.xlim(np.min(wbins), np.max(wbins))
#
## plt.axvline(x=0,color='k',linewidth=0.5)
#ax.set_xlabel(units_var2)
#
#cbar = plt.colorbar(im, ax=ax, shrink=0.75)
#cbar.ax.set_ylabel('%')
#
#plt.savefig(figdir+'cfad_vmfp_'+imemb+'_'+itest+'.png',dpi=200, facecolor='white', \
#            bbox_inches='tight', pad_inches=0.2)


## In[37]:
#
#
## create figure
#fig = plt.figure(figsize=(14,8))
#ax = fig.add_subplot(111)
#
#ax.set_title(r"$\theta_v'$")
#ax.set_ylabel('Pressure [hPa]')
#
#ax.invert_yaxis()
#ax.set_yscale('log')
#
## fill contour
#clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
#
#im = ax.contourf(thvbins[0:nbin-1], pres, np.transpose(thv_freq), clevs, locator=ticker.LogLocator(), \
#                 cmap='GnBu', alpha=1, extend='max', zorder=2)
#
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#
#plt.xlim(np.min(thvbins), np.max(thvbins))
#
#plt.axvline(x=0,color='k',linewidth=0.5)
#ax.set_xlabel('K')
#
#cbar = plt.colorbar(im, ax=ax, shrink=0.75)
#cbar.ax.set_ylabel('%')
#
#plt.show()
## plt.savefig(figdir+figtag+'_compcross_'+imemb+'_'+ivar_select+'.png',dpi=200, facecolor='white', \
##             bbox_inches='tight', pad_inches=0.2)
#
#
## In[79]:
#
#
## create figure
#fig = plt.figure(figsize=(14,8))
#ax = fig.add_subplot(111)
#
#ax.set_title(r"$\theta_v', weighted by |qv|$")
#ax.set_ylabel('Pressure [hPa]')
#
#ax.invert_yaxis()
#ax.set_yscale('log')
#
## fill contour
#clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
#
#im = ax.contourf(thvbins[0:nbin-1], pres, np.transpose(thv_freq), clevs, locator=ticker.LogLocator(), \
#                 cmap='GnBu', alpha=1, extend='max', zorder=2)
#
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#
#plt.xlim(np.min(thvbins), np.max(thvbins))
#
#plt.axvline(x=0,color='k',linewidth=0.5)
#ax.set_xlabel('K')
#
#cbar = plt.colorbar(im, ax=ax, shrink=0.75)
#cbar.ax.set_ylabel('%')
#
#plt.show()
## plt.savefig(figdir+figtag+'_compcross_'+imemb+'_'+ivar_select+'.png',dpi=200, facecolor='white', \
##             bbox_inches='tight', pad_inches=0.2)
#
#
## In[40]:
#
#
## create figure
#fig = plt.figure(figsize=(14,8))
#ax = fig.add_subplot(111)
#
#ax.set_title(r"$\theta_e'$")
#ax.set_ylabel('Pressure [hPa]')
#
#ax.invert_yaxis()
#ax.set_yscale('log')
#
## fill contour
#clevs=[0.01,0.05,0.1,0.5,1,5,10,50]
#
#im = ax.contourf(thvbins[0:nbin-1], pres, np.transpose(thv_freq), clevs, locator=ticker.LogLocator(), \
#                 cmap='GnBu', alpha=1, extend='max', zorder=2)
#
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#
#plt.xlim(np.min(thvbins), np.max(thvbins))
#
#plt.axvline(x=0,color='k',linewidth=0.5)
#ax.set_xlabel('K')
#
#cbar = plt.colorbar(im, ax=ax, shrink=0.75)
#cbar.ax.set_ylabel('%')
#
#plt.show()
## plt.savefig(figdir+figtag+'_compcross_'+imemb+'_'+ivar_select+'.png',dpi=200, facecolor='white', \
##             bbox_inches='tight', pad_inches=0.2)

