# ### Notebook to plot mean vertical profiles of vertical motion and other variables from WRF model output for new stratiform paper.
# 
# James Ruppert  
# jruppert@ou.edu  
# 7/11/24

import numpy as np
from thermo_functions import *
from precip_class import *
from memory_usage import *
from read_functions import *
import pickle
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()

print('Node ',comm.rank, 'WORKING!!')

# #### Main settings

# Testing mode: only runs a couple time steps and doesn't write out
testing=True
testing=False

# Update files instead of writing from scratch?
update_files = True
# update_files = False

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'
main_pickle = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/jruppert/tc_postproc/"+storm+'/'
figdir = "/home/jamesrup/figures/tc/ens/boxplot/"

# Members
nmem = 10 # number of ensemble members (1-5 have NCRF)
# nmem = 3
enstag = str(nmem)

# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

# Get dimensions
def get_nt(test_str):
    datdir = main+storm+'/'+memb_all[0]+'/'+test_str+'/'+datdir2
    nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
    # dp = (pres[1]-pres[0])*1e2 # Pa
    # nt=np.min([nt,nt_data-time_neglect])
    return nt_data, nz, pres
nt_ctl, nz, pres = get_nt('ctl')

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    # tests = ['ctl','ncrf36h']
    # tests = ['crfon','ncrf']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']

ntest = len(tests)

### FUNCTIONS ###############################################

# ### Functions to compute means
# These means are averaged over all ensemble members plus some selection of convective type or moisture threshold.

# def get_mean_indices(pclass, cwv):
def get_mean_indices(datdir, t0, t1):
    # Get indices for each mean
    # 0-5: pclass 0-5
    # 6: pclass=1,4,5 (deep conv + strat + anvil)
    # 7 = where PW ≥ 48 (Mapes et al. 2018)
    # 8 = where PW < 48 (Mapes et al. 2018)
    # 9 = whole domain
    mean_str = ["Non-cloud", "Deep", "Cong", "Shallow", "Strat", "Anvil",
                'MCS', 'Moist', 'Dry', 'All']
    nmean = len(mean_str)

    # Precip classification
    q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
    pclass = precip_class(q_int)
    pclass = pclass[:,np.newaxis,...] # Add single-element vertical dimension

    # CWV
    varname='PW'
    cwv = var_read_2d(datdir,varname,t0,t1,mask=True,drop=True) # mm
    cwv = cwv[:,np.newaxis,...] # Add single-element vertical dimension

    pclass_z = np.repeat(pclass, nz, axis=1)
    cwv_z    = np.repeat(cwv, nz, axis=1)

    moist_margin = 48 # kg/m2 or mm

    indices_mean_2d = {}
    indices_mean_3d = {}
    for imean in range(nmean):
        if imean <= 5:
            indices_mean_2d[mean_str[imean]] = (pclass == imean)
            indices_mean_3d[mean_str[imean]] = (pclass_z == imean)
        elif imean == 6:
            indices_mean_2d[mean_str[imean]] = ((pclass == 1) | (pclass == 4) | (pclass == 5))
            indices_mean_3d[mean_str[imean]] = ((pclass_z == 1) | (pclass_z == 4) | (pclass_z == 5))
        elif imean == 7:
            indices_mean_2d[mean_str[imean]] = (cwv >= moist_margin)
            indices_mean_3d[mean_str[imean]] = (cwv_z >= moist_margin)
        elif imean == 8:
            indices_mean_2d[mean_str[imean]] = (cwv < moist_margin)
            indices_mean_3d[mean_str[imean]] = (cwv_z < moist_margin)
        elif imean == 9:
            indices_mean_2d[mean_str[imean]] = (cwv < 1e9)
            indices_mean_3d[mean_str[imean]] = (cwv_z < 1e9)

    return mean_str, indices_mean_2d, indices_mean_3d

# ### Function to process a single test for a given ensemble member
def process_member(datdir, main_pickle, memb_str, test_str):

    # Time settings
    # Run all available time steps
    t0=0
    nt_test, nz, pres = get_nt(test_str)
    t1=nt_test

    # Read ZB once to add to Z
    zb = var_read_zb_hires(datdir, mask=True, drop=True)

    # Run only TEST time steps
    # t1_test=49
    # if test_str == 'ctl':
    #     # t0=time_neglect
    #     t0=36
    #     # Control test time sample
    #     t1=t0+t1_test
    # else:
    #     t0=0
    #     # t1=49 # max
    #     # Control test time sample
    #     t1=t1_test
    # For testing
    if testing:
        t0=30
        t1=33

    # Get indices to average over for each mean
    # mean_str, indices_mean_2d, indices_mean_3d = get_mean_indices(pclass, cwv)
    mean_str, indices_mean_2d, indices_mean_3d = get_mean_indices(datdir, t0, t1)
    nmean = len(mean_str)

    # Compute all means across all variables for this ens member

    def compute_means(mean_str, indices_mean, invar):
        nmean = len(mean_str)
        var_mean = {}
        for imean in range(nmean):
            meanvar = np.mean(invar.data, axis=(2,3), where=indices_mean[mean_str[imean]])
            var_mean[mean_str[imean]] = np.squeeze(meanvar)
        return var_mean

    def read_mean_3d_var(datdir, t0, t1, varname, mean_str, indices_mean_3d):
        var_tmp = var_read_3d_hires(datdir, varname, t0, t1, mask=True, drop=True)
        if varname == 'Z':
            var_tmp += zb
        var_mean = compute_means(mean_str, indices_mean_3d, var_tmp)
        return var_mean

    def vert_int(invar, pres, vint_top):
        # Vertically integrate: 1/g * SUM(var*dp)
        # Negative is absorbed by dp>0
        dp = pres[0]-pres[1] # Pa
        g = 9.81 # m/s2
        # vint_top = 100e2 # top for MSE integrals
        k_vint_top = np.where(pres == vint_top)[0][0]
        return np.sum(invar[:,:k_vint_top+1], axis=1)*dp/g

    def compute_vint(mean_str, invar_3d, pres, vint_top=100e2):
        nmean = len(mean_str)
        var_vint = {}
        for imean in range(nmean):
            var_vint[mean_str[imean]] = vert_int(invar_3d[mean_str[imean]], pres, vint_top)
        return var_vint

    def read_mean_vmf_vars(datdir, t0, t1, pres, mean_str, indices_mean_2d, indices_mean_3d):
        rho = var_read_3d_hires(datdir, 'rho', t0, t1, mask=True, drop=True)
        w = var_read_3d_hires(datdir, 'W', t0, t1, mask=True, drop=True)

        g = 9.81 # m/s2
        mse = var_read_3d_hires(datdir, 'mse', t0, t1, mask=True, drop=True)
        dse = var_read_3d_hires(datdir, 'dse', t0, t1, mask=True, drop=True)
        mse_vadv = -w * g*rho * np.gradient(mse, pres, axis=1)
        dse_vadv = -w * g*rho * np.gradient(dse, pres, axis=1)

        wu = np.where((w > 0), w, 0)
        wd = np.where((w < 0), w, 0)
        # wu_vint = vert_int(wu, pres, 100e2)
        # wd_vint = vert_int(wd, pres, 100e2)
        # wu_tmp = np.where(wu_vint, (wu_vint != 0), np.nan)
        # wu_tmp = np.ma.masked_where(wu_vint, (wu_vint == 0))
        # pe = 1 - (-wd_vint/wu_tmp)
        # # pe_mean = compute_means(mean_str, indices_mean_2d, pe[:,np.newaxis,...])
        # pe_mean = {}
        # nmean = len(mean_str)
        # for imean in range(nmean):
        #     mean_tmp = np.nanmean(pe[:,np.newaxis,...], axis=(2,3), where=indices_mean_2d[mean_str[imean]])
        #     pe_mean[mean_str[imean]] = np.squeeze(mean_tmp)
        #     # pe_mean[mean_str[imean]] = np.squeeze(np.ma.filled(mean_tmp, np.nan))

        # Get mean profiles
        rho_mean = compute_means(mean_str, indices_mean_3d, rho)
        w_mean = compute_means(mean_str, indices_mean_3d, w)
        wu_mean = compute_means(mean_str, indices_mean_3d, wu)
        wd_mean = compute_means(mean_str, indices_mean_3d, wd)
        mse_mean = compute_means(mean_str, indices_mean_3d, mse)
        mse_vadv_mean = compute_means(mean_str, indices_mean_3d, mse_vadv)
        dse_vadv_mean = compute_means(mean_str, indices_mean_3d, dse_vadv)
        # For doing sum over area
        # rho_mean = compute_mean_profiles(mean_str, indices_mean_3d, rho)
        # w_mean = {}
        # wu_mean = {}
        # wd_mean = {}
        # for imean in range(nmean):
        #     w_mean[mean_str[imean]]  = np.sum(w_tmp.data, axis=(2,3), where=indices_mean_3d[mean_str[imean]])
        #     wu_mean[mean_str[imean]] = np.sum(wu.data,   axis=(2,3), where=indices_mean_3d[mean_str[imean]])
        #     wd_mean[mean_str[imean]] = np.sum(wd.data,   axis=(2,3), where=indices_mean_3d[mean_str[imean]])

        # Vertical integration
        vmf = compute_vint(mean_str, w_mean, pres)
        vmfu = compute_vint(mean_str, wu_mean, pres)
        vmfd = compute_vint(mean_str, wd_mean, pres)
        # pev2 = {}
        # for imean in range(nmean):
        #     pev2[mean_str[imean]] = 1 - (vmfd[mean_str[imean]]/vmfu[mean_str[imean]])
        vmfu_600 = compute_vint(mean_str, wu_mean, pres, vint_top=600e2)
        vmfd_600 = compute_vint(mean_str, wd_mean, pres, vint_top=600e2)
        mse_vint = compute_vint(mean_str, mse_mean, pres)
        vadv_mse_vint = compute_vint(mean_str, mse_vadv_mean, pres)
        vadv_dse_vint = compute_vint(mean_str, dse_vadv_mean, pres)

        diag_vars = {
            'rho': rho_mean,
            'W': w_mean,
            'vmf': vmf,
            'vmfu': vmfu,
            'vmfd': vmfd,
            # 'pe': pe_mean,
            # 'pev2': pev2,
            'vmfu_500': vmfu_600,
            'vmfd_500': vmfd_600,
            'mse_vint': mse_vint,
            'vadv_mse_vint': vadv_mse_vint,
            'vadv_dse_vint': vadv_dse_vint,
        }

        return diag_vars

    def read_process_condh(datdir, t0, t1, pres, mean_str, indices_mean_3d):
        condh = var_read_3d_hires(datdir, 'H_DIABATIC', t0, t1, mask=True, drop=True) # K/s
        cp = 1004 # J/(kg*K)
        condh *= cp # Convert to J/(kg*s), yields W/m2 once integrated
        # Get mean profiles
        condh_mean = compute_means(mean_str, indices_mean_3d, condh)
        # Vertical integration
        # condh_vint = compute_vint(mean_str, condh_mean, pres)
        # return condh_vint
        return condh_mean

    def read_process_advec(datdir, t0, t1, mean_str, indices_mean_3d):
        dse = var_read_3d_hires(datdir, 'dse', t0, t1, mask=True, drop=True) # J/kg
        u = var_read_3d_hires(datdir, 'U', t0, t1, mask=True, drop=True) # m/s
        v = var_read_3d_hires(datdir, 'V', t0, t1, mask=True, drop=True) # m/s
        # Compute advection
        dx = 3000 # m, horizontal grid spacing
        dse_uadv = u * np.gradient(dse, axis=3)/dx # J/(kg*s)
        dse_vadv = v * np.gradient(dse, axis=2)/dx # J/(kg*s)
        # Get mean profiles
        dse_uadv_mean = compute_means(mean_str, indices_mean_3d, dse_uadv)
        dse_vadv_mean = compute_means(mean_str, indices_mean_3d, dse_vadv)
        return dse_uadv_mean, dse_vadv_mean

    def read_process_rain(datdir, t0, t1, mean_str, indices_mean_2d):
        rain = var_read_2d(datdir, 'rainrate', t0, t1, mask=True, drop=True) # mm/(time step)
        # lv0=2.5e6 # J/kg
        # rain_wm2 = rain*lv0/(24*3600) # mm/d to W/m2
        # rain_wm2 = rain*lv0/(3600) # mm/(time step) to W/m2
        # Get mean profiles
        rain_mean = compute_means(mean_str, indices_mean_2d, rain[:,np.newaxis,...])
        return rain_mean

    def read_process_condh_kik(datdir, t0, t1, mean_str, indices_mean_2d):
        varfil_main = Dataset(datdir+'condh_kukulies.nc')
        condh = varfil_main.variables['condh_k'][t0:t1,:,:]
        varfil_main.close()
        condh = mask_edges(condh, mask=True, drop=True)
        # Get mean profiles
        condh_mean = compute_means(mean_str, indices_mean_2d, condh[:,np.newaxis,...])
        return condh_mean

    varnames=[
        'QVAPOR',
        'T',
        'RTHRATLW',
        'RTHRATLWC',
        'RTHRATSW',
        'RTHRATSWC',
        'Z',
    ]

    # Pickle file to save means
    # pickle_file = main_pickle+memb_all[imemb]+'/mean_profiles_'+str(t1_test)+'hrs.pkl'
    pickle_file = main_pickle+memb_str+'/mean_profiles_'+test_str+'_alltime_v2.pkl'

    if update_files:
        # Read existing pickle file
        with open(pickle_file, 'rb') as file:
            allvars_3d_mean = pickle.load(file)
            print('Read existing pickle file: ',pickle_file)

        ##########################################################
        # Place code to process new/updated variables below here
        ##########################################################

        allvars_3d_mean['rain'] = read_process_rain(datdir, t0, t1, mean_str, indices_mean_2d)

        if testing:
            print("Test worked! Ending job before write-out...")
            sys.exit()

        ### Save updated pickle file
        with open(pickle_file, 'wb') as file:
            pickle.dump(allvars_3d_mean, file)

    else:

        # Read 3D variables and take means at same time
        allvars_3d_mean = {}

        # Special case for W/VMF
        diag_vars = read_mean_vmf_vars(datdir, t0, t1, pres*1e2, mean_str, indices_mean_2d, indices_mean_3d)
        for key in diag_vars.keys():
            allvars_3d_mean[key] = diag_vars[key]
        # for mean_name in mean_str:
        #     print('PE (',mean_name,'): ', allvars_3d_mean['pe'][mean_name])
        #     print('PEv2 (',mean_name,'): ', allvars_3d_mean['pev2'][mean_name])

        # Special case for CONDH (H_DIABATIC)
        allvars_3d_mean['condh'] = read_process_condh(datdir, t0, t1, pres*1e2, mean_str, indices_mean_3d)

        # Special case for rainfall
        allvars_3d_mean['rain'] = read_process_rain(datdir, t0, t1, mean_str, indices_mean_2d)

        # The rest of the variables
        for varname in varnames:
            allvars_3d_mean[varname] = read_mean_3d_var(datdir, t0, t1, varname, mean_str, indices_mean_3d)

        # Remove height - since don't need this
        # allvars_3d_mean.pop('Z')

        if testing:
            print("Test worked! Ending job before write-out...")
            sys.exit()

        ### Save new pickle file
        with open(pickle_file, 'wb') as file:
            pickle.dump(allvars_3d_mean, file)

    return None

# ### Read loop

# Main read loops for 3D (dependent) variables

for itest in range(ntest):
# for itest in range(1,ntest):

    test_str = tests[itest]
    print()
    print('Running test: ',test_str)
    print()

    # Loop over ensemble members
    # for imemb in range(nmem):
    # for ii in range(2):
        # imemb = comm.rank + ii*5
    imemb = comm.rank #+7

    print('Running imemb: ',memb_all[imemb])
    print()

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print('Datdir: ',datdir)
    print()

    # Process ensemble member
    process_member(datdir, main_pickle, memb_all[imemb], test_str)

print(comm.rank, 'DONE!!')
