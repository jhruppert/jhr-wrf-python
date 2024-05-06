# Script to regrid (coarsen) dataset and write out the result.
# 
# James Ruppert  
# 2/15/24

import numpy as np
# import xarray as xr
from thermo_functions import *
from precip_class import *
from read_functions import *
import xesmf as xe
from write_ncfile import *
from time import time as runtimer

################################################
#### Main settings

# Size of new grid in nx x ny
nxx = 15 # n-points in both horizontal dimensions
buffer = 80 # equal to buffer used for other routines

storm = 'haiyan'
# storm = 'maria'

filename_out='regrid_diag_nx'+str(nxx)+'.nc'

main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'

# Tests to read and compare
if storm == 'haiyan':
    # tests = ['ctl']
    # tests = ['ctl','ncrf36h']
    tests = ['crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    tests = ['ctl','ncrf36h','crfon60h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
elif storm == 'maria':
    tests = ['ctl','ncrf36h']
    tests = ['ncrf36h']
    tests = ['ctl','ncrf48h','ncrf36h']
    # tests = [tests[1],'crfon72h']
    # tests = ['crfon72h']

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 1

################################################
# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

# Get dimensions
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
nt_ctl, nz, nx1, nx2, pres = get_file_dims(datdir)
dp = (pres[1]-pres[0])*1e2 # Pa

# Get WRF file list
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
wrffiles, lat, lon = get_wrf_filelist(datdir)
lat = lat[:,0]
lon = lon[0,:]

################################################
# Output grid settings
hnxx = int(nxx/2)

# New lon/lat coords for output
x1ind = range(buffer,nx1-buffer,nxx)
x2ind = range(buffer,nx2-buffer,nxx)
nx1_new = len(x1ind)
nx2_new = len(x2ind)
lat_new = lat[x1ind]
lon_new = lon[x2ind]

# Need grid bounds for CONSERVATIVE regrid
hdy=0.5*(lat[1]-lat[0])
lat_b = np.append(lat[:]-hdy, lat[nx1-1]+hdy) # Edge points
lat_new_b = np.append(lat_new[:]-hdy, lat_new[nx1_new-1]+hdy) # Edge points
hdx=0.5*(lon[1]-lon[0])
lon_b = np.append(lon[:]-hdx, lon[nx2-1]+hdx) # Edge points
lon_new_b = np.append(lon_new[:]-hdx, lon_new[nx2_new-1]+hdx) # Edge points

# Create regridder
grid_in = {"lon": lon, "lat": lat,
           "lon_b": lon_b, "lat_b": lat_b}
grid_out = {"lon": lon_new, "lat": lat_new,
           "lon_b": lon_new_b, "lat_b": lat_new_b}
regridder = xe.Regridder(grid_in, grid_out, "conservative")

################################################
#### NetCDF variable metadata

def var_regrid_metadata(nt,nz,nx1_new,nx2_new):
    
    var_names = [
        'lat_new',
        'lon_new',
        'pres',
        'pclass_area',
        'rain',
        'qrain',
        'qrain_z',
        'qtotal',
        'pw',
        'pw_sat',
        'vmfu',
        'vmfd',
        'condh',
        'mse_vint',
        'lwacre',
        'swacre',
        'theta_e',
        'w',
        'rho',
    ]
    descriptions = [
        'latitude of new grid',
        'longitude of new grid',
        'pressure',
        'precip class area',
        'rain rate (centered diff)',
        'column integrated rain water mixr',
        'rain water mixr at lowest vertical level',
        'total integrated hydrometeor mixr',
        'precipitable water (aka CWV)',
        'saturation PW or CWV',
        'upward-masked mass flux vertically integrated (up to 100 hPa)',
        'downward-masked mass flux vertically integrated (up to 100 hPa)',
        'condensation heating from H_DIABATIC vertically int (up to 100 hPa), converted to rainfall units',
        'vertically int moist static energy, calculated as 1/g*integral(mse)dp up to 100 hPa',
        'LW column ACRE',
        'SW column ACRE',
        'equivalent potential temperature',
        'vertical motion',
        'density',
    ]
    units = [
        'deg',
        'deg',
        'hPa',
        '%',
        'mm/day',
        'mm',
        'kg/kg',
        'mm',
        'mm',
        'mm',
        'kg/m/s',
        'kg/m/s',
        'mm/day',
        'J/kg',
        'W/m^2',
        'W/m^2',
        'K',
        'm/s',
        'kg/m^3',
    ]
    dims2d = (nt,nx1_new,nx2_new)
    dims3d = (nt,nz,nx1_new,nx2_new)
    dim_names = ('nt','nx1_new','nx2_new')
    dim_names3d = ('nt','nz','nx1_new','nx2_new')
    dims_set = [
        [('nx1_new',),(nx1_new,)],
        [('nx2_new',),(nx2_new,)],
        [('nz',),(nz,)],
        [(dim_names[0],'pclass',dim_names[1],dim_names[2]), (dims2d[0],6,dims2d[1],dims2d[2])],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names,dims2d],
        [dim_names3d,dims3d],
        [dim_names3d,dims3d],
        [dim_names3d,dims3d],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set

################################################
# Function to calculate precip-class area fraction

def get_pclass_area(pclass, nt, nx1_new, nx2_new, x1ind, x2ind, nxx, hnxx):

    pclass_area = np.zeros((nt,6,nx1_new,nx2_new))

    nxsqd = np.square(nxx)
    for it in range(nt):
        for ix1 in range(nx1_new):
            for ix2 in range(nx2_new):
                x1loc=x1ind[ix1]
                x2loc=x2ind[ix2]
                for ipclass in range(6):
                    indices_pclass = (pclass[it, x1loc-hnxx:x1loc+hnxx, x2loc-hnxx:x2loc+hnxx] == ipclass).nonzero()
                    pclass_area[it,ipclass,ix1,ix2] = 100*indices_pclass[0].shape[0]/nxsqd

    return pclass_area

##### Main loops and calculations ######################################

print()
print('Running storm: ',storm)

ntest=len(tests)
for ktest in range(ntest):

    test_str=tests[ktest]

    print()
    print('Running test: ',test_str)

    # Loop over ensemble members

    for imemb in range(nmem):
    # for imemb in range(1,nmem):

        start = runtimer()

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2

        # Variable list to process

        varname='PW'
        varfil_read = Dataset(datdir+varname+'.nc')
        pw = np.squeeze(varfil_read.variables[varname][:,:,:,:])
        varfil_read.close()
        nt=pw.shape[0]

        t0=0
        t1=nt

        var_names, descriptions, units, dims_set = var_regrid_metadata(nt,nz,nx1_new,nx2_new)

        # Stratiform ID
        q_int = read_qcloud(datdir,t0,t1,mask=False) # mm
        pclass = precip_class(q_int)

        # MSE diag vars
        # PW Saturated
        pw_sat = read_mse_diag(datdir,'pw_sat',2,t0,t1,mask=False)
        # VMFU, VMFD
        vmfu = read_mse_diag(datdir,'vmfu',2,t0,t1,mask=False)
        vmfd = read_mse_diag(datdir,'vmfd',2,t0,t1,mask=False)
        # MSE (vertically integrated)
        mse = read_mse_diag(datdir,'mse_vint',2,t0,t1,mask=False)
        # Microphysics heating
        condh = read_mse_diag(datdir,'condh',2,t0,t1,mask=False)

        # Rain rate
        varname = 'rainrate'
        rain = var_read_2d(datdir,varname,t0,t1,mask=False) # mm/d

        # QINT (including QRAIN) column integrated
        qcloud  = q_int[0] # mm
        qrain  = q_int[1] # mm
        qice  = q_int[2] # mm
        qsnow  = q_int[3] # mm
        qgraup  = q_int[4] # mm
        qtotal = qcloud + qrain + qice + qsnow + qgraup
        
        # QRAIN: lowest level
        varfil_main = Dataset(datdir+'QRAIN'+'.nc')
        qrain_lowest = np.squeeze(varfil_main.variables['QRAIN'][t0:t1,0,:,:])
        varfil_main.close()

        # ACRE
        lwacre = read_lwacre(datdir,t0,t1,mask=False) # W/m2
        swacre = read_swacre(datdir,t0,t1,mask=False) # W/m2

        # Equiv potential temp
        varname = 'QVAPOR'
        qv = var_read_3d(datdir,varname,t0,t1,mask=False) # kg/kg
        varname = 'T'
        tmpk = var_read_3d(datdir,varname,t0,t1,mask=False) # K
        theta_e = theta_equiv(tmpk,qv,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

        # Density
        rho = density_moist(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # kg/m3

        # Vertical motion
        varname = 'W'
        w = var_read_3d(datdir,varname,t0,t1,mask=False) # m/s

        # ### Interpolate variable ##############################################

        print("Running regridding...")

        var_list=[]
        var_list.append(lat_new)
        var_list.append(lon_new)
        var_list.append(pres)

        print("  ... running pclass area")
        pclass_area = get_pclass_area(pclass, nt, nx1_new, nx2_new, x1ind, x2ind, nxx, hnxx)
        var_list.append(pclass_area)
        print("  ... running regridder")
        var_list.append(regridder(rain.data))
        var_list.append(regridder(qrain.data))
        var_list.append(regridder(qrain_lowest.data))
        var_list.append(regridder(qtotal.data))
        var_list.append(regridder(pw.data))
        var_list.append(regridder(pw_sat.data))
        var_list.append(regridder(vmfu.data))
        var_list.append(regridder(vmfd.data))
        var_list.append(regridder(condh.data))
        var_list.append(regridder(mse.data))
        var_list.append(regridder(lwacre.data))
        var_list.append(regridder(swacre.data))
        var_list.append(regridder(theta_e.data))
        var_list.append(regridder(w.data))
        var_list.append(regridder(rho.data))

        # ### Write out variables ##############################################

        print("Writing netcdf")
        write_ncfile(datdir+filename_out, var_list, var_names, descriptions, units, dims_set)

        end = runtimer()
        time_elapsed = end - start
        print("Time elapsed for member: ", time_elapsed)
