
#  ### Notebook to genereate cross sections from TC output binned according to a 2D variable.
# 
#  Assumes output is in a single netcdf file on pressure levels.
# 
#  James Ruppert
#  jruppert@ou.edu
#  4/23/22


# NOTE: Using copied tracking from CTL for NCRF tests

import numpy as np
from thermo_functions import *
from precip_class import precip_class
from memory_usage import memory_usage
from read_functions import *
import pickle
import sys
# from mpi4py import MPI

# comm = MPI.COMM_WORLD
# nproc = comm.Get_size()


#  #### Main settings


# Write out pickle file?
# do_write=True
# do_write=False

do_hires=True # Vertical high resolution?
# do_hires=False

# do_tests=True # Run sensitivity tests?
do_tests=False


# Index variable (2D; independent var)
# ivar_select = 'pw'
ivar_select = 'sf'
# ivar_select = 'lwacre'
# ivar_select = 'rain'
# options (requiring 2D info): pw, rain, lwacre
# options (requiring 3D info): vmf

# Fill variable (3D; dependent var)
# fillvar_select = 'lwcrf'
# fillvar_select = 'h_diabatic'
# options: avor, lwcrf, tprm, dbz, rh

# Contour variable (3D; dependent var)
contvar_select = 'w'

# Mask out all points except [stratiform/nonrain/etc], or switch off
# istrat=-1 #2 # 0-non-raining, 1-conv, 2-strat, 3-other/anvil, (-1 for off)
# istrat=2

# Number of sample time steps
nt=200 # will be chopped down to max available
# nt=48
# nt=6
time_neglect=12 # time steps from start to neglect (for CTL only)
t1_test=12 # n time steps to sample for tests
# t1_test=6
t1_test=24
# t1_test=48 # didn't work -- too much memory


#  #### Additional settings and directories


# storm = 'haiyan'
storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
main_pickle = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/jruppert/tc_postproc/"+storm+'/'
datdir2 = 'post/d02/'

# Time selection
# hr_tag = str(np.char.zfill(str(nt), 2))

# Tests to read and compare
# tests = ['crfon','ncrf']
if storm == 'haiyan':
    if do_tests:
        tests = ['ctl','ncrf36h','STRATANVIL_OFF','STRATANVIL_ON','STRAT_OFF']
    else:
        tests = ['ctl']
    # tests = ['ctl']
    t0_test=36
elif storm == 'maria':
    if do_tests:
        # tests = ['ctl','ncrf36h']
        tests = ['ctl','ncrf48h']
    else:
        tests = ['ctl']
    t0_test=48
ntest=len(tests)

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 2
enstag = str(nmem)

# Shift starting-read time step for CRFON comparison
# t0_test=0
# if tests[0] == 'crfon':
#     t0_test=24
#     memb0=5 # for CRFFON test

# TC tracking
# ptrack='600' # tracking pressure level
# var_track = 'rvor' # variable
# rmax = 8 # radius (deg) limit for masking around TC center

# Strat/Conv index subset
# if istrat == -1:
#     fig_extra=''
#     strattag=''
# else:
#     if istrat == 0:
#         strattag='Nonrain'
#     elif istrat == 1:
#         strattag='Conv'
#     elif istrat == 2:
#         strattag='Strat'
#     elif istrat == 3:
#         strattag='Anv'
#     fig_extra='_'+strattag.lower()



# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

# Get dimensions
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'+datdir2
nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
dp = (pres[1]-pres[0])*1e2 # Pa
# if do_tests:
nt=t1_test
nt=np.min([nt,nt_data-time_neglect])

# For dropped edges
nx1-=80*2
nx2-=80*2
if do_hires:
    nz=39
    pres = np.arange(1000,25,-25)

# Get WRF file list
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
wrffiles, lat, lon = get_wrf_filelist(datdir)



#  #### Main loops and compositing


# Main read loops for 3D (dependent) variables

for ktest in range(ntest):
    # ktest = comm.rank
    test_str=tests[ktest]
    print('Running test: ',test_str)

    # if do_write:

    # Arrays to save variables
    dims = (nmem,nt,nz,nx1,nx2)
    dims2d = (nmem,nt,nx1,nx2)
    crf_all = np.ma.zeros(dims, dtype=np.float32)
    hdia_all = np.ma.zeros(dims, dtype=np.float32)
    tmpk_all = np.ma.zeros(dims, dtype=np.float32)
    qv_all = np.ma.zeros(dims, dtype=np.float32)
    cvar_all = np.ma.zeros(dims, dtype=np.float32)
    strat_all = np.ma.zeros(dims2d, dtype=np.float32)
    rain_all = np.ma.zeros(dims2d, dtype=np.float32)
    # ddtq_all = np.ma.zeros(dims2d, dtype=np.float32)
    # mse_all = np.ma.zeros(dims2d)
    # mselw_all = np.ma.zeros(dims2d)
    lwacre_all = np.ma.zeros(dims2d, dtype=np.float32)
    sh_all = np.ma.zeros(dims2d, dtype=np.float32)
    lh_all = np.ma.zeros(dims2d, dtype=np.float32)

    ivar_all = np.ma.zeros(dims2d)

    # t0=time_neglect # neglect the first 12 time steps
    # t1=t0+nt
    if test_str == 'ctl':
        t0=time_neglect
        t1=nt+t0
        if do_tests:
            t0=t0_test
            # t1=t0+49
            # Control test time sample
            t1=t0+t1_test
    else:
        t0=0
        # t1=49 # max
        # Control test time sample
        t1=t1_test

    # t0+=1 # add one time step since NCRF(t=0) = CTL
    # t1 = t0+nt

    # Loop over ensemble members

    for imemb in range(nmem):

        print('Running imemb: ',memb_all[imemb])

        datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
        datdir3d = datdir #+'v2/'
        print(datdir)

        # Localize to TC track
        # NOTE: Using copied tracking from CTL for NCRF tests
        # track_file = datdir+'../../track_'+var_track+'_'+ptrack+'hPa.nc'
        # trackfil_ex=''
        # if 'crf' in test_str:
        #     trackfil_ex='_ctlcopy'
        # track_file = datdir+'../../track_'+var_track+trackfil_ex+'_'+ptrack+'hPa.nc'

        # Required variables

        # Stratiform index
        q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
        strat_all[imemb,:,:,:] = precip_class(q_int)
        del q_int

        # MSE variance
        # mse = read_mse(datdir,t0,t1) # J/m2
        lwacre_all[imemb,:,:,:] = read_lwacre(datdir,t0,t1,mask=True,drop=True) # W/m2

        # Surface fluxes
        varname = 'HFX'
        sh_all[imemb,:,:,:] = var_read_2d(datdir,varname,t0,t1,mask=True, drop=True) # W/m2
        varname = 'LH'
        lh_all[imemb,:,:,:] = var_read_2d(datdir,varname,t0,t1,mask=True, drop=True) # W/m2

        # Line contour variable ("cvar")
        # Vertical motion
        if contvar_select == 'w':
            varname='W'
            # cvar_all[ktest,imemb,:,:,:,:] = var_read_3d(datdir3d,varname,t0,t1)*1e2 # m/s --> cm/s
            cvar_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True)*1e2 # m/s --> cm/s
            units_var2='cm/s'
            lcmin = -20; lcmax=20; lcint=2

        # CWV
        varname='PW'
        # cwv = var_read_2d(datdir3d,varname,t0,t1) # mm
        # ddtq = np.gradient(cwv, axis=0) # mm/hr

        # Rain rate
        varname = 'rainrate'
        rain_all[imemb,...] = var_read_2d(datdir,varname,t0,t1,mask=True,drop=True) # mm/d

        # Three-dimensional dependent variables ("var")

        varname = 'RTHRATLWCRF'
        crf_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) * 3600.*24 # K/s --> K/d
        varname = 'H_DIABATIC'
        hdia_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) * 3600.*24 # K/s --> K/d
        varname = 'T'
        tmpk_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) # K
        varname = 'QVAPOR'
        qv_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) # kg/kg
        # # Radar Reflectivity
        # if fillvar_select == 'dbz':
        #     varname = fillvar_select
        #     var = var_read_3d(datdir3d,varname,t0,t1)
        # # Radiation
        # elif fillvar_select == 'lwcrf':
        #     varname = 'RTHRATLWCRF'
        #     var_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) * 3600.*24 # K/s --> K/d
        #     # varname = 'RTHRATLWC'
        #     # var -= var_read_3d(datdir3d,varname,t0,t1) * 3600.*24 # K/s --> K/d
        # # Horizontal temperature anomaly
        # elif fillvar_select == 'tprm':
        #     varname = 'T'
        #     tmpk = var_read_3d(datdir3d,varname,t0,t1) # K
        #     varname = 'QVAPOR'
        #     qv = var_read_3d(datdir3d,varname,t0,t1) # kg/kg
        #     var = theta_virtual(tmpk,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K
        #     # Subtract time-dependent domain average
        #     radius_ls = 12 # deg
        #     # var_ls = mask_tc_track(track_file, radius_ls, var, lon, lat, t0, t1)
        #     # var_ls_avg = np.ma.mean(var_ls,axis=(0,2,3))
        #     # var -= var_ls_avg[np.newaxis,:,np.newaxis,np.newaxis]
        # # Relative humidity
        # elif fillvar_select == 'rh':
        #     varname = 'T'
        #     tmpk = var_read_3d(datdir3d,varname,t0,t1) # K
        #     # print(np.min(tmpk))
        #     varname = 'QVAPOR'
        #     qv = var_read_3d(datdir3d,varname,t0,t1) # kg/kg
        #     # var = relh(qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2,tmpk,ice=1) # RH in %
        # # Absolute vorticity
        # elif fillvar_select == 'avor':
        #     varname = 'AVOR'
        #     var = var_read_3d(datdir3d,varname,t0,t1) # 10^-5 /s
        #     var*=10 # --> 10^/6 /s
        # # Diabatic heating
        # elif fillvar_select == 'h_diabatic':
        #     varname = 'H_DIABATIC'
        #     var_all[imemb,:,:,:,:] = var_read_3d_hires(datdir3d,varname,t0,t1,mask=True,drop=True) * 3600.*24 # K/s --> K/d


        # Indexing variable

        # PW
        if ivar_select == 'pw':
            varname = ivar_select.upper()
            ivar_all[imemb,:,:,:] = np.squeeze(var_read_2d(datdir,varname,t0,t1,mask=True,drop=True))
        # Saturation fraction
        elif ivar_select == 'sf':
            ipw = var_read_2d(datdir,'PW',t0,t1,mask=True,drop=True)
            ipw_sat = read_mse_diag(datdir,'pw_sat',t0,t1,mask=True,drop=True)
            # ivar_all[ktest,imemb,:,:,:] = np.squeeze(100*ipw/ipw_sat)
            ivar_all[imemb,:,:,:] = np.squeeze(ipw/ipw_sat)
        # Rainfall rate
        elif ivar_select == 'rain':
            varname = 'rainrate'
            # rain = var_read_2d(datdir,varname,t0,t1) # mm/hr
            irain = var_read_3d_ik(datdir,'QRAIN',t0,t1,ik=0,mask=True,drop=True)
            ivar_all[imemb,:,:,:] = irain # kg/kg
        # LW-ACRE
        elif ivar_select == 'lwacre':
            ivar_all[imemb,:,:,:] = read_lwacre(datdir,t0,t1,mask=True,drop=True) # W/m2
        # Vertical mass flux
        elif ivar_select == 'vmf':
            g=9.81 # gravity, m/s2
            varname='W'
            w = var_read_3d(datdir3d,varname,t0,t1) # m/s
            wv_int = np.sum(w,axis=1) * dp/g # m/s * s**2/m * kg/m/s**2 = kg/s/m
            ivar_all[imemb,:,:,:] = np.reshape(wv_int,(nt,1,nx1,nx2))
        # Theta-e (equivalent potential temperature)
        # elif ivar_select == 'th_e':
        #     varname='T'
        #     tmpk = var_read_3d(datdir3d,varname,t0,t1) # K
        #     varname = 'QVAPOR'
        #     qv = var_read_3d(datdir3d,varname,t0,t1) # kg/kg
        #     th_e = theta_equiv(tmpk,qv,qv,(pres[np.newaxis,:,np.newaxis,np.newaxis])*1e2) # K

        ### Process and save variable ##############################################

        # Mask out based on strat/conv
        # if istrat != -1:
        #     var = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), var, copy=True)
        #     cvar = np.ma.masked_where((np.repeat(strat,nz,axis=1) != istrat), cvar, copy=True)
        #     sys.exit()

        # Calculate MSE LW variance term
        # rmax_mean = 10
        # mse_mn = mask_tc_track(track_file, rmax_mean, mse[:,np.newaxis,:,:], lon, lat, t0, t1)
        # lwn_mn = mask_tc_track(track_file, rmax_mean, lwnet[:,np.newaxis,:,:], lon, lat, t0, t1)
        # mse_mn = np.ma.mean(mse_mn,axis=(1,2,3))
        # lwn_mn = np.ma.mean(lwn_mn,axis=(1,2,3))
        # msep = np.squeeze(mask_tc_track(track_file, rmax, mse[:,np.newaxis,:,:], lon, lat, t0, t1))
        # msep -= mse_mn[:,np.newaxis,np.newaxis]
        # lwnp = np.squeeze(mask_tc_track(track_file, rmax, lwnet[:,np.newaxis,:,:], lon, lat, t0, t1))
        # lwnp -= lwn_mn[:,np.newaxis,np.newaxis]
        # mselw = msep * lwnp

        # Localize to TC track
        # var = mask_tc_track(track_file, rmax, var, lon, lat, t0, t1)
        # cvar = mask_tc_track(track_file, rmax, cvar, lon, lat, t0, t1)
        # strat = mask_tc_track(track_file, rmax, strat, lon, lat, t0, t1)
        # lwacre = mask_tc_track(track_file, rmax, lwacre, lon, lat, t0, t1)
        # mse = mask_tc_track(track_file, rmax, mse[:,np.newaxis,:,:], lon, lat, t0, t1)

        # Save ens member
        # var_all[ktest,imemb,:,:,:,:]  = var
        # cvar_all[ktest,imemb,:,:,:,:] = cvar
        # strat_all[ktest,imemb,:,:,:]  = strat[:,0,:,:]
        # ddtq_all[ktest,imemb,:,:,:]   = ddtq[:,0,:,:]
        # mse_all[ktest,imemb,:,:,:]    = mse[:,0,:,:]
        # mselw_all[ktest,imemb,:,:,:]  = mselw
        # lwacre_all[ktest,imemb,:,:,:] = lwacre[:,0,:,:]

    #### Calculate basic mean
    # var_mn = np.ma.mean(var_all,axis=(1,2,4,5))

    #### Calculate basic mean
    # mselw_mn = np.ma.mean(mselw_all,axis=(1,2,3,4))
    # msevar1 = np.var(mse_all, axis=(3,4)) # (J/m2)^2 --> (ktest,nmemb,nt)
    # msevar = np.mean(msevar1, axis=(1,2)) # (J/m2)^2 --> (ktest)



    # memory_usage()



    #  #### Index aka Bin variable


    # Variable settings

    # PW
    if ivar_select == 'pw':
        fmin=35;fmax=80 # mm
        step=1
        bins=np.arange(fmin,fmax+step,step)
        xlabel='Column water vapor [mm]'
        log_x='linear'
    # Column saturation fraction
    elif ivar_select == 'sf':
        # fmin=30;fmax=100 # %
        fmin=.3;fmax=1 # %
        step=0.01
        bins=np.arange(fmin,fmax+step,step)
        # xlabel='Column saturation fraction [%]'
        xlabel='Column saturation fraction'
        log_x='linear'
    # Rainfall rate
    elif ivar_select == 'rain':
        # bins=10.**(np.arange(1,8,0.3)-4)
        bins=10.**(np.arange(0,8,0.3)-4)
        xlabel='Rainfall rate [mm/hr]'
        log_x='log'
    # LW-ACRE
    elif ivar_select == 'lwacre':
        # fmin=-50; fmax=200 # W/m2
        fmin=50; fmax=200 # W/m2
        step=5
        bins=np.arange(fmin,fmax+step,step)
        xlabel='LW-ACRE [W/m**2]'
        log_x='linear'
    # Vertical mass flux
    elif ivar_select == 'vmf':
        bins=10.**(np.arange(1,8,0.3)-3)
        # bins=np.flip(-1.*bins)
        xlabel='Vertical mass flux [kg/m/s]'
        log_x='log'
    # Theta-e (equivalent potential temperature)
    elif ivar_select == 'th_e':
        fmin=315; fmax=365 # K
        step=1
        bins=np.arange(fmin,fmax+step,step)
        xlabel=r'$\theta_e$ [K]'
        log_x='linear'

    nbins = np.size(bins)

    # Create axis of bin center-points for plotting
    bin_axis = (bins[np.arange(nbins-1)]+bins[np.arange(nbins-1)+1])/2


    #  #### Conduct compositing

    # Loop and composite variables

    frequency = np.zeros(nbins-1)
    # var_binned=np.ma.zeros((nbins-1,nz))
    crf_binned=np.ma.zeros((nbins-1,nz))
    hdia_binned=np.ma.zeros((nbins-1,nz))
    tmpk_binned=np.ma.zeros((nbins-1,nz))
    qv_binned=np.ma.zeros((nbins-1,nz))
    cvar_binned=np.ma.zeros((nbins-1,nz))
    strat_binned=np.ma.zeros((nbins-1,6)) # Bin count: 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    # ddtq_binned=np.ma.zeros((ntest,nbins-1))
    # mselw_binned=np.ma.zeros((ntest,nbins-1))
    lwacre_binned=np.ma.zeros((nbins-1))
    sh_binned=np.ma.zeros((nbins-1))
    lh_binned=np.ma.zeros((nbins-1))
    rain_binned=np.ma.zeros((nbins-1))
    # mselw_strat=np.ma.zeros((ntest,nbins-1,4)) # LSEPLWP in: 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    lwacre_strat=np.ma.zeros((nbins-1,6)) # LWACRE in: 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    rain_strat=np.ma.zeros((nbins-1,6)) # rain in: 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    # ddtq_strat=np.ma.zeros((ntest,nbins-1,6)) # DDTQ in: 0-non-raining, 1-conv, 2-strat, 3-other/anvil

    # Bin the variables, averaging across member, time, x, y: (ntest,nmemb,nt,nz,nx1,nx2) --> (ntest,nbins,nz)
    for ibin in range(nbins-1):

        indices = ((ivar_all >= bins[ibin]) & (ivar_all < bins[ibin+1])).nonzero()
        frequency[ibin] = indices[0].shape[0]

        if indices[0].shape[0] > 3:
            # var_binned[ibin,:]  = np.mean(var_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            crf_binned[ibin,:]  = np.mean(crf_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            hdia_binned[ibin,:]  = np.mean(hdia_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            tmpk_binned[ibin,:]  = np.mean(tmpk_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            qv_binned[ibin,:]  = np.mean(qv_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            cvar_binned[ibin,:] = np.mean(cvar_all[indices[0],indices[1],:,indices[2],indices[3]], axis=0, dtype=np.float64)
            # ddtq_binned[ktest,ibin] = np.mean(ddtq_all[ ktest,indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
            # mselw_binned[ktest,ibin]  = np.mean(mselw_all[ktest,indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
            lwacre_binned[ibin]  = np.mean(lwacre_all[indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
            sh_binned[ibin]  = np.mean(sh_all[indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
            lh_binned[ibin]  = np.mean(lh_all[indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
            rain_binned[ibin]  = np.mean(rain_all[indices[0],indices[1],  indices[2],indices[3]], axis=0, dtype=np.float64)
        else:
            # var_binned[ibin,:] = np.nan
            crf_binned[ibin,:] = np.nan
            hdia_binned[ibin,:] = np.nan
            tmpk_binned[ibin,:] = np.nan
            qv_binned[ibin,:] = np.nan
            cvar_binned[ibin,:] = np.nan
            # ddtq_binned[ktest,ibin] = np.nan
            # mselw_binned[ktest,ibin] = np.nan
            lwacre_binned[ibin]= np.nan
            sh_binned[ibin]= np.nan
            lh_binned[ibin]= np.nan
            rain_binned[ibin] = np.nan
            strat_binned[ibin,:] = np.nan
            # mselw_strat[ktest,ibin,:] = np.nan
            lwacre_strat[ibin,:] = np.nan
            rain_strat[ibin,:] = np.nan
            # ddtq_strat[ktest,ibin,:] = np.nan
            continue

        for kstrat in range(0,6):
            indices_strat = ((ivar_all >= bins[ibin]) & (ivar_all < bins[ibin+1]) & 
            (strat_all == kstrat)).nonzero()
            # indices_strat = (strat_all[ktest,indices[0],indices[1],indices[2],indices[3]] == kstrat).nonzero()

            strat_binned[ibin,kstrat] = indices_strat[0].shape[0]

            # Bin the 2D var by rain class
            if indices_strat[0].shape[0] > 3:
                # mselw_strat[ktest,ibin,kstrat] = np.mean(mselw_all[ktest,indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0, dtype=np.float64)
                lwacre_strat[ibin,kstrat] = np.mean(lwacre_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0, dtype=np.float64)
                rain_strat[ibin,kstrat] = np.mean(rain_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0, dtype=np.float64)
                # lwacre_strat[ktest,ibin,kstrat] = np.mean(lwacre_all[ktest].flatten()[indices_strat[0]], axis=0, dtype=np.float64)
                # ddtq_strat[ktest,ibin,kstrat] = np.mean(ddtq_all[ktest,indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0, dtype=np.float64)
            else:
        #         mselw_strat[ktest,ibin,kstrat] = np.nan
                lwacre_strat[ibin,kstrat] = np.nan
                rain_strat[ibin,kstrat] = np.nan
                # ddtq_strat[ktest,ibin,kstrat] = np.nan



    #  ### Write out to pickle file


    # if do_write:
    # for itest in range(ntest):
    pickle_file = main_pickle+'/binned_2d_'+test_str+'_'+ivar_select+'_'+str(nmem)+'memb_'+str(nt)+'hrs.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump([bins, frequency, crf_binned, hdia_binned, tmpk_binned, qv_binned,
                     cvar_binned, lwacre_binned, sh_binned, lh_binned, rain_binned,
                        strat_binned, lwacre_strat, rain_strat], file)
    print('Done writing '+test_str+'!')


