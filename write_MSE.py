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
import sys
# from thermo_functions import theta_equiv, density_moist


# #### Main settings

pres_top = 100

storm = 'haiyan'
storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h']
    # tests = [tests[1],'crfon60h']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    # tests = [tests[1],'crfon72h']

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 2
enstag = str(nmem)
# Starting member to read
memb0=1


nums=np.arange(memb0,nmem+memb0,1); nums=nums.astype(str)
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

# Get Lat/Lon
process = subprocess.Popen(['ls '+main+storm+'/'+memb_all[0]+'/'+tests[0]+'/wrfout_d02_*'],shell=True,
    stdout=subprocess.PIPE,universal_newlines=True)
output = process.stdout.readline()
wrffil = output.strip() #[3]
varfil_main = Dataset(wrffil)
lat = varfil_main.variables['XLAT'][:][0] # deg
lon = varfil_main.variables['XLONG'][:][0] # deg
varfil_main.close()
lon1d=lon[0,:]
lat1d=lat[:,0]

ikread = np.where(pres == pres_top)[0][0]


# #### NetCDF variable read functions

def var_read(datdir,varname,ikread):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][:,0:ikread+1,:,:]
    varfil_main.close()
    return var


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


##### MSE / DSE convergence functions ####################

def calc_vadv(w, rho, var, dp, g):
    # Gradient terms (Inoue and Back 2015):
    #   VADV_SUM = < omeg * dVAR/dp >
    omeg = w * (-1)*g*rho
    vadv = omeg * np.gradient(var,(-dp),axis=1)
    vadv_sum = np.sum(vadv, axis=1)*dp/g
    return vadv_sum

def mse_hflux(u, v, x1d, y1d, dse, mse, dp, g):
    # Gradient terms (Inoue and Back 2015):
    #   dse-term = del . <sV> where s = DSE
    #       = d/dx <su> + d/dy <sv>
    #   mse-term = same but with h = MSE
    #   < > is vertical integral over the troposphere
    su = np.sum(u * dse, axis=1)*dp/g
    sv = np.sum(v * dse, axis=1)*dp/g
    hu = np.sum(u * mse, axis=1)*dp/g
    hv = np.sum(v * mse, axis=1)*dp/g
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
    grad_s = np.sum(dse * div, axis=1)*dp/g
    grad_h = np.sum(mse * div, axis=1)*dp/g
    return grad_s, grad_h


# #### Main loops and calculations

# Read base-state height once per storm
varname='ZB'
z_b = var_read(datdir,varname,ikread) # m

# Main read loops for 3D (dependent) variables

ntest=2

for ktest in range(ntest):

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
        cons = dp/g
        mse_int = np.sum(mse, axis=1) * cons # J/m^2

        # Diagnostics for Gross Moist Stability

        varfil_main = Dataset(datdir+'density.nc')
        rho = varfil_main.variables['rho'][:,0:ikread+1,:,:] # kg/m3
        varfil_main.close()
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

        write_vars(datdir,nt,nz,nx1,nx2, dse,mse,mse_int,
                   grad_s_vadv,   grad_h_vadv,
                   grad_s_hflux,  grad_h_hflux,
                   grad_s_converg,grad_h_converg)