#!/usr/bin/env python
# coding: utf-8

# ### Python script to write out moist static energy and its
#       vertical integral from TC output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/24/23


# NOTE: Using copied tracking from CTL for NCRF tests

from netCDF4 import Dataset
import numpy as np
import subprocess
import os
from thermo_functions import density_moist
import sys
# from thermo_functions import theta_equiv, density_moist


# #### Main settings

pres_top = 100 # top for MSE integrals

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
    # tests = ['STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    # tests = [tests[1],'crfon72h']

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
nmem = 1

######################################################################

# Ens member strings
memb0=1 # Starting member to read
nums=np.arange(memb0,nmem+memb0,1)
nums=nums.astype(str)
nustr = np.char.zfill(nums, 2)
memb_all=np.char.add('memb_',nustr)

datdir2 = 'post/d02/'

##### Get dimensions

datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
datdir3d = datdir #+'v2/'
varfil_main = Dataset(datdir3d+'T.nc')
nz = varfil_main.dimensions['level'].size
nx1 = varfil_main.dimensions['lat'].size
nx2 = varfil_main.dimensions['lon'].size
pres = varfil_main.variables['pres'][:] # hPa
dp = (pres[0]-pres[1])*1e2 # Pa
varfil_main.close()

ikread = np.where(pres == pres_top)[0][0]

# Get Lat/Lon
wrfdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
dirlist = os.listdir(wrfdir)
subs="wrfout_d02"
wrf_files = list(filter(lambda x: subs in x, dirlist))
wrf_files = [wrfdir + s for s in wrf_files]
varfil_main = Dataset(wrf_files[0])
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()
lon1d=lon[0,:]
lat1d=lat[:,0]

######################################################################

# #### NetCDF variable read functions

def var_read(datdir,varname,ikread):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,0:ikread+1,:,:]
    varfil_main.close()
    return var

######################################################################

# #### NetCDF variable write function

def write_vars(datdir,nt,nz,nx1,nx2,dse,mse,mse_int,
                grad_s_vadv,   grad_h_vadv,
                grad_s_hflux,  grad_h_hflux,
                grad_s_converg,grad_h_converg):

    file_out = datdir+'mse.nc'
    ncfile = Dataset(file_out,mode='w', clobber=True)

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    dse_nc = ncfile.createVariable('dse', np.single, ('nt','nz','nx1','nx2',))
    dse_nc.units = 'J/kg'
    dse_nc.long_name = 'dry static energy, calculated as cpT + gz'
    dse_nc[:,:,:,:] = dse[:,:,:,:]

    mse_nc = ncfile.createVariable('mse', np.single, ('nt','nz','nx1','nx2',))
    mse_nc.units = 'J/kg'
    mse_nc.long_name = 'moist static energy, calculated as cpT + gz + L_v*q'
    mse_nc[:,:,:,:] = mse[:,:,:,:]

    msei_nc = ncfile.createVariable('mse_int', np.single, ('nt','nx1','nx2',))
    msei_nc.units = 'J/m^2'
    msei_nc.long_name = 'integrated moist static energy, calculated as 1/g*integral(mse)dp up to 100 hPa'
    msei_nc[:,:,:] = mse_int[:,:,:]

    grad_s_v = ncfile.createVariable('grad_s_vadv', np.single, ('nt','nx1','nx2',))
    grad_s_v.units = 'J/m^2/s'
    grad_s_v.long_name = 'integrated VADV of DSE (up to 100 hPa)'
    grad_s_v[:,:,:] = grad_s_vadv[:,:,:]

    grad_h_v = ncfile.createVariable('grad_h_vadv', np.single, ('nt','nx1','nx2',))
    grad_h_v.units = 'J/m^2/s'
    grad_h_v.long_name = 'integrated VADV of MSE (up to 100 hPa)'
    grad_h_v[:,:,:] = grad_h_vadv[:,:,:]

    grad_s_h = ncfile.createVariable('grad_s_hflux', np.single, ('nt','nx1','nx2',))
    grad_s_h.units = 'J/m^2/s'
    grad_s_h.long_name = 'integrated del.(DSE*V) (up to 100 hPa)'
    grad_s_h[:,:,:] = grad_s_hflux[:,:,:]

    grad_h_h = ncfile.createVariable('grad_h_hflux', np.single, ('nt','nx1','nx2',))
    grad_h_h.units = 'J/m^2/s'
    grad_h_h.long_name = 'integrated del.(MSE*V) (up to 100 hPa)'
    grad_h_h[:,:,:] = grad_h_hflux[:,:,:]

    grad_s_c = ncfile.createVariable('grad_s_converg', np.single, ('nt','nx1','nx2',))
    grad_s_c.units = 'J/m^2/s'
    grad_s_c.long_name = 'integrated DSE(del.V) (up to 100 hPa)'
    grad_s_c[:,:,:] = grad_s_converg[:,:,:]

    grad_h_c = ncfile.createVariable('grad_h_converg', np.single, ('nt','nx1','nx2',))
    grad_h_c.units = 'J/m^2/s'
    grad_h_c.long_name = 'integrated MSE(del.V) (up to 100 hPa)'
    grad_h_c[:,:,:] = grad_h_converg[:,:,:]

    ncfile.close()


def var_metadata():
    
    var_names = [
        # 'dse',
        'mse',
        # 'mse_vint',
        # 'grad_s_vadv',
        # 'grad_h_vadv',
        'grad_s_hflux',
        'grad_h_hflux',
        'grad_s_converg',
        'grad_h_converg',
    ]
    long_names = [
        # 'dry static energy, calculated as cpT + gz',
        'moist static energy, calculated as cpT + gz + L_v*q',
        # 'integrated moist static energy, calculated as 1/g*integral(mse)dp up to 100 hPa',
        # 'integrated VADV of DSE (up to 100 hPa)',
        # 'integrated VADV of MSE (up to 100 hPa)',
        'integrated del.(DSE*V) (up to 100 hPa)',
        'integrated del.(MSE*V) (up to 100 hPa)',
        'integrated DSE(del.V) (up to 100 hPa)',
        'integrated MSE(del.V) (up to 100 hPa)',
    ]

    units = [
        'J/kg',
        # 'J/kg',
        # 'J/m^2',
        # 'J/m^2/s',
        # 'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
        'J/m^2/s',
    ]

    len1 = len(var_names); len2 = len(long_names); len3 = len(units)
    if (len1 != len2) or (len1 != len3):
        raise ValueError("Variable info counts are off")

    return var_names, long_names, units

def new_write_vars(file_out, var_list, var_names, long_names, units):

    ncfile = Dataset(file_out,mode='w', clobber=True)

    nt, nz, nx1, nx2 = var_list[0].shape

    time_dim = ncfile.createDimension('nt', nt) # unlimited axis (can be appended to).
    z_dim = ncfile.createDimension('nz', nz)
    x1_dim = ncfile.createDimension('nx1', nx1)
    x2_dim = ncfile.createDimension('nx2', nx2)

    dims4d = ('nt','nz','nx1','nx2')
    dims3d = ('nt','nx1','nx2')

    nvar = len(var_list)
    for ivar in range(nvar):
        vardims = len(var_list[ivar].shape)
        if vardims == 3:
            dimsout = dims3d
        else:
            dimsout = dims4d
        writevar = ncfile.createVariable(var_names[ivar], np.single, dimsout)
        writevar.units = units[ivar]
        writevar.long_name = long_names[ivar]
        writevar[...] = var_list[ivar]

    ncfile.close()

##### MSE / DSE convergence functions ####################

def vert_int(var, dp, g):
    # Vertically integrate: 1/g * SUM(var*dp)
    # Assumes dp is a constant
    # Negative is absorbed by dp>0
    var_sum = np.sum(var, axis=1)*dp/g
    return var_sum

def calc_vadv(w, rho, var, dp, g):
    # Gradient terms (Inoue and Back 2015):
    #   VADV_SUM = < omeg * dVAR/dp >
    omeg = w * (-1)*g*rho
    vadv = omeg * np.gradient(var,(-dp),axis=1)
    vadv_sum = vert_int(vadv, dp, g)
    return vadv_sum

def mse_hflux(u, v, x1d, y1d, dse, mse, dp, g):
    # Gradient terms (Inoue and Back 2015):
    #   dse-term = del . <sV> where s = DSE
    #       = d/dx <su> + d/dy <sv>
    #   mse-term = same but with h = MSE
    #   < > is vertical integral over the troposphere
    su = vert_int((u * dse), dp, g)
    sv = vert_int((v * dse), dp, g)
    hu = vert_int((u * mse), dp, g)
    hv = vert_int((v * mse), dp, g)
    grad_s_x = np.gradient(su,x1d,axis=2)
    grad_s_y = np.gradient(sv,y1d,axis=1)
    grad_s = grad_s_x + grad_s_y
    grad_h_x = np.gradient(hu,x1d,axis=2)
    grad_h_y = np.gradient(hv,y1d,axis=1)
    grad_h = grad_h_x + grad_h_y
    return grad_s, grad_h

def mse_converg(u, v, x1d, y1d, dse, mse, dp, g):
    # Gradient terms (Inoue and Back 2015):
    #   dse-term = <s del . V> where s = DSE
    #   mse-term = same but with h = MSE
    dudx = np.gradient(u,x1d,axis=3) # /s
    dvdy = np.gradient(v,y1d,axis=2) # /s
    div = dudx + dvdy
    grad_s = vert_int((dse * div), dp, g)
    grad_h = vert_int((mse * div), dp, g)
    return grad_s, grad_h


# #### Main loops and calculations

# Read base-state height once per storm
varname='ZB'
z_b = var_read(datdir,varname,ikread) # m

# Main read loops for 3D (dependent) variables

ntest=len(tests)
# for ktest in range(ntest):
for ktest in range(1,2):

    test_str=tests[ktest]

    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    
        print('Running imemb: ',memb_all[imemb])
    
        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        print(datdir)

        # Required variables

        # Water vapor
        varname='QVAPOR'
        qv = var_read(datdir,varname,ikread) # kg/kg
        nt,nz,nx1,nx2 = qv.shape
        # Temperature
        varname='T'
        tmpk = var_read(datdir,varname,ikread) # K
        # Height
        varname='Z'
        z = var_read(datdir,varname,ikread) # m2/s2
        z += z_b

        # Calculate Moist Static Energy (MSE)

        # ;LATENT HEAT OF VAPORIZATION
        cp=1004.  # J/K/kg
        cpl=4186. # J/k/kg
        cpv=1885. # J/K/kg
        lv0=2.5e6 # J/kg
        lv = lv0 - (cpl-cpv)*(tmpk-273.15)

        g=9.81 # m/s2
        dse = cp*tmpk + g*z
        mse = dse + lv*qv

        # Vertically integrate
        mse_int = vert_int(mse, dp, g) # J/m^2

        # Diagnostics for Gross Moist Stability
        rho = density_moist(tmpk,qv,pres[np.newaxis,0:ikread+1,np.newaxis,np.newaxis,]*1e2) # kg/m3
        w = var_read(datdir,'W',ikread) # m/s
        grad_s_vadv = calc_vadv(w, rho, dse, dp, g)
        grad_h_vadv = calc_vadv(w, rho, mse, dp, g)

        deg2m = np.pi*6371*1e3/180
        x1d = lon1d * deg2m
        y1d = lat1d * deg2m
        u = var_read(datdir,'U',ikread)
        v = var_read(datdir,'V',ikread)

        grad_s_hflux, grad_h_hflux     = mse_hflux(  u, v, x1d, y1d, dse, mse, dp, g)
        grad_s_converg, grad_h_converg = mse_converg(u, v, x1d, y1d, dse, mse, dp, g)

        ### Write out variables ##############################################

        var_list=[]
        var_list.append(mse)
        var_list.append(grad_s_hflux)
        var_list.append(grad_h_hflux)
        var_list.append(grad_s_converg)
        var_list.append(grad_h_converg)

        var_names, long_names, units = var_metadata()

        file_out = datdir+'testout.nc'
        new_write_vars(file_out, var_list, var_names, long_names, units)
        # write_vars(datdir,nt,nz,nx1,nx2, dse,mse,mse_int,
        #            grad_s_vadv,   grad_h_vadv,
        #            grad_s_hflux,  grad_h_hflux,
        #            grad_s_converg,grad_h_converg)