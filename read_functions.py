# Functions to read in and process variables

from netCDF4 import Dataset
import numpy as np
import subprocess

###### Functions for getting file specs #########################

def get_file_dims(datdir):
    # varfil_main = Dataset(datdir+'T_HiRes.nc')
    varfil_main = Dataset(datdir+'AVOR_HiRes.nc')
    # varfil_main = Dataset(datdir+'U10.nc')
    nz = varfil_main.dimensions['level'].size
    nx1 = varfil_main.dimensions['lat'].size
    nx2 = varfil_main.dimensions['lon'].size
    nt = varfil_main.dimensions['time'].size
    pres = varfil_main.variables['pres'][:] # hPa
    varfil_main.close()
    dims = (nt, nz, nx1, nx2)
    return nt, nz, nx1, nx2, pres

def get_wrf_filelist(datdir):
    process = subprocess.Popen(['ls '+datdir+'wrfout_d02_*'],shell=True,
        stdout=subprocess.PIPE,universal_newlines=True)
    wrffiles = process.stdout.readlines()
    for ifil in range(len(wrffiles)):
        wrffiles[ifil] = wrffiles[ifil].strip()
    varfil_main = Dataset(wrffiles[0])
    lat = varfil_main.variables['XLAT'][:][0] # deg
    lon = varfil_main.variables['XLONG'][:][0] # deg
    varfil_main.close()
    return wrffiles, lat, lon


###### Function to mask edges #########################

def mask_edges(array, mask=True, drop=False):
    # 
    # Will now just return unmodified array if mask=False
    # 
    # Returns a masked array with edges masked
    # Last dimensions of input array must be x1,x2
    #   It is otherwise versatile
    buffer=80
    # DROP = True will lop off edges instead of masking them
    if mask:
        if drop:
            nx1 = array.shape[-2]
            nx2 = array.shape[-1]
            array = array[...,buffer:nx2-buffer]
            array = array[...,buffer:nx1-buffer,:]
        else:
            array = np.ma.array(array, mask=False, copy=False)
            array[...,0:buffer,:]=np.ma.masked
            array[...,-buffer:,:]=np.ma.masked
            array[...,:,0:buffer]=np.ma.masked
            array[...,:,-buffer:]=np.ma.masked
        # array = np.ma.filled(array, fill_value=np.nan)
    return array


###### Functions for reading variables #########################

def var_read_3d(datdir, varname, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return mask_edges(var,mask,drop)

def var_read_3d_ik(datdir, varname, t0, t1, ik, mask=True, drop=False):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][t0:t1,ik,:,:]
    varfil_main.close()
    return mask_edges(var,mask,drop)

def var_read_3d_hires(datdir, varname, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+varname+'_HiRes.nc')
    var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return mask_edges(var,mask,drop)

def var_read_2d(datdir, varname, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+varname+'.nc')
    var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return mask_edges(np.squeeze(var),mask,drop)


###### Special variable reads ###################################

def read_rain(datdir, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+'rainrate_HiRes.nc')
    rain = varfil_main.variables['RAINC'][t0:t1,:,:,:]
    varfil_main.close()
    return mask_edges(np.squeeze(rain),mask,drop)

def read_qcloud(datdir, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+'q_int.nc')
    q_int = varfil_main.variables['q_int'][:,t0:t1,:,:]
    varfil_main.close()
    return mask_edges(np.squeeze(q_int),mask,drop)

# def read_mse(datdir,t0,t1,drop=False):
#     varfil_main = Dataset(datdir+'mse.nc')
#     mse = varfil_main.variables['mse_int'][t0:t1,:,:] # J/m2
#     varfil_main.close()
#     return mask_edges(mse,drop)

def var_read_zb_hires(datdir, mask=True, drop=False):
    varfil_main = Dataset(datdir+'ZB_HiRes.nc')
    var = varfil_main.variables['ZB'][:,:,:,:]
    varfil_main.close()
    return mask_edges(var,mask,drop)

def read_mse_diag(datdir, varname, t0, t1, mask=True, drop=False):
    varfil_main = Dataset(datdir+'mse_diag.nc')
    dims = varfil_main.variables[varname].dimensions
    ndims = len(dims)
    if ndims == 3: # (t,y,x)
        var = varfil_main.variables[varname][t0:t1,:,:]
    elif ndims == 4: # (t,p,y,x)
        var = varfil_main.variables[varname][t0:t1,:,:,:]
    varfil_main.close()
    return mask_edges(np.squeeze(var),mask,drop)

def read_lwnet(datdir, t0, t1, mask=True, drop=False):
    lw_t = var_read_2d(datdir,'LWUPT',t0,t1,mask,drop) - var_read_2d(datdir,'LWDNT',t0,t1,mask,drop) # W/m2
    lw_b = var_read_2d(datdir,'LWUPB',t0,t1,mask,drop) - var_read_2d(datdir,'LWDNB',t0,t1,mask,drop) # W/m2
    lwnet = lw_b - lw_t
    return lwnet

def read_lwnetc(datdir, t0, t1, mask=True, drop=False):
    lw_t = var_read_2d(datdir,'LWUPTC',t0,t1,mask,drop) - var_read_2d(datdir,'LWDNTC',t0,t1,mask,drop) # W/m2
    lw_b = var_read_2d(datdir,'LWUPBC',t0,t1,mask,drop) - var_read_2d(datdir,'LWDNBC',t0,t1,mask,drop) # W/m2
    lwnetc = lw_b - lw_t
    return lwnetc

def read_lwacre(datdir, t0, t1, mask=True, drop=False):
    lwnet = read_lwnet(datdir,t0,t1,mask,drop) # W/m2
    lwnetc = read_lwnetc(datdir,t0,t1,mask,drop) # W/m2
    lwacre = lwnet - lwnetc
    return lwacre

def read_swnet(datdir, t0, t1, mask=True, drop=False):
    sw_t = var_read_2d(datdir,'SWUPT',t0,t1,mask,drop) - var_read_2d(datdir,'SWDNT',t0,t1,mask,drop) # W/m2
    sw_b = var_read_2d(datdir,'SWUPB',t0,t1,mask,drop) - var_read_2d(datdir,'SWDNB',t0,t1,mask,drop) # W/m2
    swnet = sw_b - sw_t
    return swnet

def read_swnetc(datdir, t0, t1, mask=True, drop=False):
    sw_t = var_read_2d(datdir,'SWUPTC',t0,t1,mask,drop) - var_read_2d(datdir,'SWDNTC',t0,t1,mask,drop) # W/m2
    sw_b = var_read_2d(datdir,'SWUPBC',t0,t1,mask,drop) - var_read_2d(datdir,'SWDNBC',t0,t1,mask,drop) # W/m2
    swnetc = sw_b - sw_t
    return swnetc

def read_swacre(datdir, t0, t1, mask=True, drop=False):
    swnet = read_swnet(datdir,t0,t1,mask,drop) # W/m2
    swnetc = read_swnetc(datdir,t0,t1,mask,drop) # W/m2
    swacre = swnet - swnetc
    return swacre
