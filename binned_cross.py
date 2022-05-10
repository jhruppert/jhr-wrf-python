#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
#import subprocess

#### Time selection

t0 = 48
t1 = t0+24

storm = 'haiyan'


#### Directories

# main = "/Users/jamesruppert/code/tc_output/"
#main = "/Users/jruppert/code/tc_output/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
figdir = "/home/jamesrup/idl/figures/tc/ens"
#storm = subprocess.check_output("ls $main", shell=True, capture_output=True)
# print(storm)
#istorm=storm[0]

# Ensemble members
#memb = !ls $main/$istorm
nums=np.arange(1,21,1); nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)
imemb=memb_all[0]
# print(main+istorm+'/'+imemb)

datdir = main+storm+'/'+imemb+'/ctl/'
datdir = datdir+'post/d02/'
print(datdir)


#### Read variables

# Fill contour variable

# Radar Reflectivity
# varfil_main = Dataset(datdir+'dbz.nc') # this opens the netcdf file
# binvar_f_in = varfil_main.variables['dbz'][t0:t1,:,:,:]
# title = 'Ref'
# units_var1 = 'dBZ'
# cmin = -20; cmax=20

# Radiation
varfil_main = Dataset(datdir+'RTHRATLW.nc') # this opens the netcdf file
binvar_f_in = varfil_main.variables['RTHRATLW'][t0:t1,:,:,:] * 3600.*24 # K/s --> K/d
varcs = Dataset(datdir+'RTHRATLWC.nc') # this opens the netcdf file
cs = varcs.variables['RTHRATLWC'][t0:t1,:,:,:] * 3600.*24 # K/s --> K/d
binvar_f_in -= cs
title = 'LW-CRF'
units_var1 = 'K/d'
cmax=4; cmin=-1.*cmax

# Conv/strat separation: varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain
varfil_strat = Dataset(datdir+'strat.nc') # this opens the netcdf file
strat_in = varfil_strat.variables['strat'][t0:t1,:,:,:]

bv_shape = np.shape(binvar_f_in)
print("Binvar shape: ",bv_shape)
nt = bv_shape[0]
nz = bv_shape[1]


# Vertical coordinate
pres = varfil_main.variables['pres'][:] # Pa
print("Vertical shape: ",np.shape(pres))


# Line contour variable

# Vertical motion
varfil_cvar = Dataset(datdir+'W.nc') # this opens the netcdf file
binvar_c_in = varfil_cvar.variables['W'][t0:t1,:,:,:]*1e2 # m/s --> mm/s
units_var2='cm/s'
lcmin = -20; lcmax=20; lcint=2


# Indexing variable
binfil = Dataset(datdir+'PW.nc') # this opens the netcdf file
ivar = binfil.variables['PW'][t0:t1,:,:,:]
print("PW shape: ",np.shape(ivar))


#### Create PW Bins

fmin=35 # mm
fmax=80
step=1
bins = np.arange(fmin,fmax,step)
nbins = np.size(bins)

#### Bin the target variable

binvar_f = np.zeros((nbins,nt,nz)) # nbins, nt, nz
binvar_c = np.zeros((nbins,nt,nz))
binvar_strat = np.zeros((nbins,nt,3))

# for ibin in range(nbins):
for itim in range(nt):
    for ibin in range(nbins):
        indices = ((ivar[itim,0,:,:] >= bins[ibin]-0.5*step) & (ivar[itim,0,:,:] < bins[ibin]+0.5*step)).nonzero()
        tmp_f = binvar_f_in[itim,:,indices[0],indices[1]]
        binvar_f[ibin,itim,:] = np.mean(tmp_f,axis=0,dtype=np.float64)
        tmp_c = binvar_c_in[itim,:,indices[0],indices[1]]
        binvar_c[ibin,itim,:] = np.mean(tmp_c,axis=0,dtype=np.float64)
        tmp_strat = strat_in[itim,:,indices[0],indices[1]]
        for istrat in range(3):
            iindex = ((tmp_strat == (istrat+1))).nonzero()
            # print(iindex)
            binvar_strat[ibin,itim,0] = np.shape(iindex)[1]

#### Time-average

binvar_f_mn = np.mean(binvar_f,axis=1)
binvar_c_mn = np.mean(binvar_c,axis=1)
binvar_s_mn = np.mean(binvar_strat,axis=1)



### Plotting routines

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


### CONTOUR PLOT

# create figure
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111)

ax.set_title('Binned plot')
ax.set_ylabel('Pressure [hPa]')
ax.set_xlabel('Column water vapor [mm]')

# fill contour
nlevs=21
inc=(cmax-cmin)/nlevs
clevs = np.arange(cmin, cmax+inc, inc)
pltvar=binvar_f_mn
im = ax.contourf(bins, pres, np.transpose(pltvar), clevs, cmap='RdBu_r', alpha=0.6, \
                 extend='max', zorder=2)

cbar = plt.colorbar(im, ax=ax, shrink=0.75)
cbar.ax.set_ylabel(units_var1)
ax.invert_yaxis()
ax.set_yscale('log')
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())

# ax2=ax.twinx()
# im = ax.plot(bins, binvar_s_mn)

# line contour
clevs = [1,2,5,10,50,100,500,1000,2000,3000]#np.arange(lcmin, lcmax, lcint)
clevs = np.concatenate((-1*np.flip(clevs),clevs))
cpltvar=binvar_c_mn
# cpltvar=np.gradient(binvar_c_mn,10000,axis=1,edge_order=2)*-1e5
im = ax.contour(bins, pres, np.transpose(cpltvar), clevs, colors='black', zorder=2)
ax.clabel(im, im.levels, inline=True, fontsize=13)

# plt.show()
plt.savefig(figdir+'myplot.png',dpi=200)


### LINE PLOT

# create figure
fig = plt.figure(figsize=(14,4))
ax = fig.add_subplot(111)

ax.set_title('Binned plot')
ax.set_ylabel('Cell ID')
ax.set_xlabel('Column water vapor [mm]')

# Conv/strat separation: varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain

pltvar=binvar_s_mn

plt.plot(bins, binvar_s_mn[:,0], "-b", label="sine")
plt.plot(bins, binvar_s_mn[:,1], "-r", label="cosine")
plt.plot(bins, binvar_s_mn[:,2], "-g", label="cosine")

plt.legend(loc="upper left")
# plt.ylim(-1.5, 2.0)

#plt.show()
plt.savefig(figdir+'myplot2.png',dpi=200)

