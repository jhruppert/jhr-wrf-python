# ### Notebook to genereate plots of binned 2D variables.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/22/24

# NOTE: Using copied tracking from CTL for NCRF tests

import numpy as np
import matplotlib.pyplot as plt
from precip_class import precip_class
from read_functions import *
import pickle

# #### Main settings

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
figdir = "/home/jamesrup/figures/tc/ens/"+storm+'/'
datdir2 = 'post/d02/'
save_dir = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/james/binned_2d_sav/"

# Time selection
# hr_tag = str(np.char.zfill(str(nt), 2))

# Tests to read and compare
# tests = ['crfon','ncrf']
# if storm == 'haiyan':
#     tests = ['ctl','ncrf36h']
# elif storm == 'maria':
#     # tests = ['ctl','ncrf36h']
#     tests = ['ctl','ncrf48h']
tests = ['ctl']

# Number of sample time steps
nt=200 # will be chopped down to max available
nt=3*24#24#6
time_neglect=12 # time steps from start to neglect

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 2
enstag = str(nmem)


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
nt=np.min([nt,nt_data-time_neglect])
nx1-=80*2
nx2-=80*2

# Get WRF file list
datdir = main+storm+'/'+memb_all[0]+'/'+tests[0]+'/'
wrffiles, lat, lon = get_wrf_filelist(datdir)


# Main read loops for all variables

# Arrays to save variables
# ntest=len(tests)
dims2d = (nmem,nt,nx1,nx2)
pclass_all = np.zeros(dims2d)
tqi_all = np.zeros(dims2d)
pw_all = np.zeros(dims2d)
satfrac_all = np.zeros(dims2d)
lwacre_all = np.zeros(dims2d)
rain_all = np.zeros(dims2d)
# vmfd_all = np.zeros(dims2d)

t0=time_neglect # neglect the first 12 time steps
t1=t0+nt

ktest=0
test_str=tests[ktest]
print('Running test: ',test_str)

# Loop over ensemble members

for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print(datdir)

    # Precip class
    q_int = read_qcloud(datdir,t0,t1,drop=True) # mm
    pclass_all[imemb,:,:,:] = precip_class(q_int)

    # IWP
    tqi_all[imemb,:,:,:] = q_int[2] # Just cloud ice
    # tqi_all[imemb,:,:,:] = q_int[2] + q_int[3] + q_int[4] # Ice water path = ice + snow + graupel

    # LWACRE
    lwacre_all[imemb,:,:,:] = read_lwacre(datdir,t0,t1,drop=True) # W/m2

    # Rain rate
    varname = 'rainrate'
    rain_all[imemb,:,:,:] = var_read_2d(datdir,varname,t0,t1,drop=True)/24 # mm/d --> mm/hr

    # PW, Sat frac
    # pw = var_read_2d(datdir3d,varname,t0,t1,drop=True) # mm
    pw = read_mse_diag(datdir,'pw',2,t0,t1,drop=True) # mm
    pw_sat = read_mse_diag(datdir,'pw_sat',2,t0,t1,drop=True) # mm
    pw_all[imemb,:,:,:] = pw
    satfrac_all[imemb,:,:,:] = 100*pw/pw_sat

    # VMFD
    # vmfd_all[imemb,:,:,:] = read_mse_diag(datdir,'vmfd',2,t0,t1,drop=True) # kg/m/s


# Bin variable settings

def binvar_settings(ivar_select, pw_all, satfrac_all, rain_all, lwacre_all):

    nbins=30

    # PW
    if ivar_select == 'pw':
        ivar_all = pw_all
        fmin=35;fmax=80 # mm
        # step=1
        bins=np.linspace(fmin,fmax,num=nbins)
        xlabel='Column Water Vapor [mm]'
        log_x='linear'
    # Column saturation fraction
    elif ivar_select == 'sf':
        ivar_all = satfrac_all
        fmin=30;fmax=102 # %
        # step=2
        bins=np.linspace(fmin,fmax,num=nbins)
        xlabel='Saturation Fraction [%]'
        log_x='linear'
    # Rainfall rate
    elif ivar_select == 'rain':
        ivar_all = rain_all
        # bins=10.**(np.arange(1,8,0.3)-4)
        # bins=10.**(np.arange(0,8,0.3)-4)
        bins=np.logspace(-4,2.5,num=nbins)
        xlabel='Rainfall Rate [mm/hr]'
        log_x='log'
    # LW-ACRE
    elif ivar_select == 'lwacre':
        ivar_all = lwacre_all
        fmin=-50; fmax=200 # W/m2
        # step=5
        bins=np.linspace(fmin,fmax,num=nbins)
        xlabel='LW-ACRE [W/m**2]'
        log_x='linear'
    # Stratiform area fraction
    # elif ivar_select == 'strat_area':
    #     fmin=0;fmax=60 # %
    #     step=1
    #     bins=np.arange(fmin,fmax+step,step)
    #     xlabel='Stratiform area fraction [%]'
    #     log_x='linear'

    # Create axis of bin center-points for plotting
    # nbins = np.size(bins)
    bin_axis = (bins[np.arange(nbins-1)]+bins[np.arange(nbins-1)+1])/2

    return ivar_all, bins, bin_axis, xlabel, log_x


# Binning function

def run_dual_binning(bins_x, bins_y, ivar_x, ivar_y, pclass_all, pw_all, satfrac_all, lwacre_all, rain_all, tqi_all):

    # Loop and composite variables

    nbins_x = np.size(bins_x)
    nbins_y = np.size(bins_y)

    nclass=6

    bin_freq=np.zeros((nbins_x-1, nbins_y-1)) # Bin counts

    pclass_binned=np.full((nbins_x-1,nbins_y-1,nclass), np.nan) # Bin count: 0-non-raining, 1-conv, 2-strat, 3-other/anvil
    pw_binned=np.full((nbins_x-1,nbins_y-1), np.nan)
    satfrac_binned=np.full((nbins_x-1,nbins_y-1), np.nan)
    lwacre_binned=np.full((nbins_x-1,nbins_y-1), np.nan)
    rain_binned=np.full((nbins_x-1,nbins_y-1), np.nan)
    tqi_binned=np.full((nbins_x-1,nbins_y-1), np.nan)
    # vmfd_binned=np.full((nbins-1), np.nan)

    pw_class=np.full((nbins_x-1,nbins_y-1,nclass), np.nan) # Binned by precip_class: 0-non-raining, 1-deep conv, 2-congest, 3-shallow, 4-strat, 5-anvil
    satfrac_class=np.full((nbins_x-1,nbins_y-1,nclass), np.nan)
    lwacre_class=np.full((nbins_x-1,nbins_y-1,nclass), np.nan)
    rain_class=np.full((nbins_x-1,nbins_y-1,nclass), np.nan)
    tqi_class=np.full((nbins_x-1,nbins_y-1,nclass), np.nan)
    # vmfd_class=np.full((nbins-1,nclass), np.nan)

    nmin = 3

    # Bin the variables, averaging across member, time, x, y: (ntest,nmemb,nt,nz,nx1,nx2) --> (ntest,nbins,nz)
    for ibin_x in range(nbins_x-1):
        for ibin_y in range(nbins_y-1):

            indices = ((ivar_x >= bins_x[ibin_x]) & (ivar_x < bins_x[ibin_x+1]) & 
                       (ivar_y >= bins_y[ibin_y]) & (ivar_y < bins_y[ibin_y+1])).nonzero()

            ifreq = indices[0].shape[0]
            bin_freq[ibin_x,ibin_y] = ifreq

            if ifreq > nmin:
                pw_binned[ibin_x,ibin_y]       = np.mean(pw_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
                satfrac_binned[ibin_x,ibin_y]  = np.mean(satfrac_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
                lwacre_binned[ibin_x,ibin_y]   = np.mean(lwacre_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
                rain_binned[ibin_x,ibin_y]     = np.mean(rain_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
                tqi_binned[ibin_x,ibin_y]      = np.mean(tqi_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
                # vmfd_binned[ibin]     = np.mean(vmfd_all[indices[0],indices[1],indices[2],indices[3]], axis=0)
            else:
                continue
                # Else will leave bins filled with NaN

            for kclass in range(nclass):
                indices_strat = ((ivar_x >= bins_x[ibin_x]) & (ivar_x < bins_x[ibin_x+1]) &
                                 (ivar_y >= bins_y[ibin_y]) & (ivar_y < bins_y[ibin_y+1]) &
                                 (pclass_all == kclass)).nonzero()

                ifreq = indices_strat[0].shape[0]
                pclass_binned[ibin_x,ibin_y,kclass] = ifreq

                # Bin the 2D var by rain class
                if ifreq > nmin:
                    pw_class[ibin_x,ibin_y,kclass] = np.mean(pw_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)
                    satfrac_class[ibin_x,ibin_y,kclass] = np.mean(satfrac_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)
                    lwacre_class[ibin_x,ibin_y,kclass] = np.mean(lwacre_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)
                    rain_class[ibin_x,ibin_y,kclass] = np.mean(rain_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)
                    tqi_class[ibin_x,ibin_y,kclass] = np.mean(tqi_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)
                    # vmfd_class[ibin,kclass] = np.mean(vmfd_all[indices_strat[0],indices_strat[1],indices_strat[2],indices_strat[3]], axis=0)

    # Calculate Area Fraction for each class
    pclass_area=np.ma.zeros((nbins_x-1,nbins_y-1,nclass))
    total=np.nansum(pclass_binned)
    for kclass in range(nclass):
        pclass_area[:,:,kclass] = pclass_binned[:,:,kclass]/total*1e2

    binned_vars = {
        'bins_x':bins_x, 'bins_y':bins_y,
        'bin_freq':bin_freq, 'pclass_binned':pclass_binned, 'pclass_area':pclass_area, 'pw_binned':pw_binned, 'satfrac_binned':satfrac_binned,
        'lwacre_binned':lwacre_binned, 'rain_binned':rain_binned, 'tqi_binned':tqi_binned, #'vmfd_binned':vmfd_binned,
        'pw_class':pw_class, 'satfrac_class':satfrac_class, 'lwacre_class':lwacre_class, 'rain_class':rain_class, 'tqi_class':tqi_class}#,
        # 'vmfd_class':vmfd_class, }

    return binned_vars


# ### Run binning and plotting

def write_pickle(save_file, binned_vars):
    with open(save_file, 'wb') as f:
        pickle.dump(binned_vars, f)


# ivar_select_x='rain'
# ivar_select_y='sf'

# ivar_x, bins_x, bin_axis_x, xlabel, log_x = binvar_settings(ivar_select_x, pw_all, satfrac_all, rain_all, lwacre_all)
# ivar_y, bins_y, bin_axis_y, ylabel, log_y = binvar_settings(ivar_select_y, pw_all, satfrac_all, rain_all, lwacre_all)

# binned_vars = run_dual_binning(bins_x, bins_y, ivar_x, ivar_y, pclass_all, pw_all, satfrac_all, lwacre_all, rain_all, tqi_all)

# save_file=save_dir+ivar_select_x+'-'+ivar_select_y+'.pkl'
# write_pickle(save_file, binned_vars)


# ivar_select_x='sf'
# ivar_select_y='lwacre'

# ivar_x, bins_x, bin_axis_x, xlabel, log_x = binvar_settings(ivar_select_x, pw_all, satfrac_all, rain_all, lwacre_all)
# ivar_y, bins_y, bin_axis_y, ylabel, log_y = binvar_settings(ivar_select_y, pw_all, satfrac_all, rain_all, lwacre_all)

# binned_vars = run_dual_binning(bins_x, bins_y, ivar_x, ivar_y, pclass_all, pw_all, satfrac_all, lwacre_all, rain_all, tqi_all)

# save_file=save_dir+ivar_select_x+'-'+ivar_select_y+'.pkl'
# write_pickle(save_file, binned_vars)


ivar_select_x='rain'
ivar_select_y='lwacre'

ivar_x, bins_x, bin_axis_x, xlabel, log_x = binvar_settings(ivar_select_x, pw_all, satfrac_all, rain_all, lwacre_all)
ivar_y, bins_y, bin_axis_y, ylabel, log_y = binvar_settings(ivar_select_y, pw_all, satfrac_all, rain_all, lwacre_all)

binned_vars = run_dual_binning(bins_x, bins_y, ivar_x, ivar_y, pclass_all, pw_all, satfrac_all, lwacre_all, rain_all, tqi_all)

save_file=save_dir+ivar_select_x+'-'+ivar_select_y+'.pkl'
write_pickle(save_file, binned_vars)
