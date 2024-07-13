# ### Notebook to plot mean vertical profiles of vertical motion and other variables from WRF model output for new stratiform paper.
# 
# James Ruppert  
# jruppert@ou.edu  
# 7/11/24

import numpy as np
import matplotlib
from matplotlib import ticker
import matplotlib.pyplot as plt
import sys
from thermo_functions import *
from precip_class import *
from memory_usage import *
from read_functions import *
import pickle
import os.path
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()


# #### Main settings

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'
main_pickle = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/jruppert/tc_postproc/"+storm+'/'
figdir = "/home/jamesrup/figures/tc/ens/boxplot/"

# Set to true to do sensitivty test comparisons
# Else: analysis of CTL only
do_tests=True
# do_tests=False
t1_test=24 # n time steps to sample for tests
# t1_test=2 # n time steps to sample for tests

# Read and write variables, or load in saved pickle?
# do_write=True
# do_write=False

# pickle_file = main_pickle+'mean_profiles_'+str(t1_test)+'hrs.pkl'
# if not do_write:
#     if not os.path.isfile(pickle_file):
#         raise Exception('Pickle file not found!')

time_neglect=12 # time steps from start to neglect

# Number of sample time steps (if only running CTL)
# nt=200 # will be chopped down to max available
nt=24
# nt=12

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
test_str='ctl'
datdir = main+storm+'/'+memb_all[0]+'/'+test_str+'/'+datdir2
nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
dp = (pres[1]-pres[0])*1e2 # Pa
nt=np.min([nt,nt_data-time_neglect])
if do_tests:
    nt=t1_test
nx1-=80*2
nx2-=80*2
# Setting for (vertical) HiRes output
nz=39
print(comm.rank, 'WORKED')

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','ncrf36h','STRATANVIL_ON','STRATANVIL_OFF','STRAT_OFF']
    tests_str = ['CTL','NCRF','STRATANVON','STRATANVOFF','STRATOFF']
    # tests = ['ctl','ncrf36h']
    # tests_str = ['CTL','NCRF']
    # tests = ['crfon','ncrf']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    tests_str = ['CTL','NCRF']

ntest = len(tests)

# ### Function to compute means
# 
# These means are averaged over all ensemble members plus some selection of convective type or moisture threshold.

def compute_mean_profiles(nmean, pclass, cwv, w, rho, qv, tmpk, lw, lwc, sw, swc):

    # 0-5: pclass 0-5
    # 6: pclass=1,4,5 (deep conv + strat + anvil)
    # 7 = where PW ≥ 48 (Mapes et al. 2018)
    # 8 = where PW < 48 (Mapes et al. 2018)
    # 9 = whole domain

    moist_margin = 48 # kg/m2 or mm

    pclass_z = np.repeat(pclass[:,np.newaxis,...], nz, axis=1)
    cwv_z    = np.repeat(cwv[:,np.newaxis,...], nz, axis=1)

    dims_mean = (nmean, nt, nz)
    w_mean    = np.zeros(dims_mean)
    rho_mean  = np.zeros(dims_mean)
    qv_mean   = np.zeros(dims_mean)
    tmpk_mean = np.zeros(dims_mean)
    lw_mean   = np.zeros(dims_mean)
    lwc_mean  = np.zeros(dims_mean)
    sw_mean   = np.zeros(dims_mean)
    swc_mean  = np.zeros(dims_mean)

    for imean in range(nmean):

        if imean <= 5:
            indices_mean = (pclass_z == imean)
        elif imean == 6:
            indices_mean = ((pclass_z == 1) | (pclass_z == 4) | (pclass_z == 5))
        elif imean == 7:
            indices_mean = (cwv_z >= moist_margin)
        elif imean == 8:
            indices_mean = (cwv_z < moist_margin)
        elif imean == 9:
            indices_mean = (cwv_z < 1e9)

        w_mean[imean,:,:]    = np.mean(w.data, axis=(2,3), where=indices_mean)
        rho_mean[imean,:,:]  = np.mean(rho.data, axis=(2,3), where=indices_mean)
        qv_mean[imean,:,:]   = np.mean(qv.data, axis=(2,3), where=indices_mean)
        tmpk_mean[imean,:,:] = np.mean(tmpk.data, axis=(2,3), where=indices_mean)
        lw_mean[imean,:,:]   = np.mean(lw.data, axis=(2,3), where=indices_mean)
        lwc_mean[imean,:,:]  = np.mean(lwc.data, axis=(2,3), where=indices_mean)
        sw_mean[imean,:,:]   = np.mean(sw.data, axis=(2,3), where=indices_mean)
        swc_mean[imean,:,:]  = np.mean(swc.data, axis=(2,3), where=indices_mean)

    return w_mean, rho_mean, qv_mean, tmpk_mean, lw_mean, lwc_mean, sw_mean, swc_mean

# ### Read loop

# Main read loops for 3D (dependent) variables

# Arrays to save variables
nmean = 10
dims = (ntest, nmean, nt, nz)
w_mean    = np.zeros(dims)
rho_mean  = np.zeros(dims)
qv_mean   = np.zeros(dims)
tmpk_mean = np.zeros(dims)
lw_mean   = np.zeros(dims)
lwc_mean  = np.zeros(dims)
sw_mean   = np.zeros(dims)
swc_mean  = np.zeros(dims)

for itest in range(ntest):
# for itest in range(1):

    test_str = tests[itest]
    print()
    print('Running test: ',test_str)

    # t0=time_neglect # neglect the first 12 time steps
    # t1=t0+nt
    if test_str == 'ctl':
        t0=time_neglect
        t1=nt+t0
        if do_tests:
            t0=36
            # t1=t0+49
            # Control test time sample
            t1=t0+t1_test
    else:
        t0=0
        # t1=49 # max
        # Control test time sample
        t1=t1_test

    # Loop over ensemble members
    # for imemb in range(nmem):
    # for imemb in range(1):
    imemb = comm.rank

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print(datdir)

    # Stratiform ID
    q_int = read_qcloud(datdir,t0,t1,mask=True,drop=True) # mm
    pclass = precip_class(q_int)

    # CWV
    varname='PW'
    cwv = var_read_2d(datdir,varname,t0,t1,mask=True,drop=True) # mm

    # 3D variables
    w = var_read_3d_hires(datdir, 'W', t0, t1, mask=True, drop=True) # m/s
    rho = var_read_3d_hires(datdir, 'rho', t0, t1, mask=True, drop=True) # kg/m3
    qv = var_read_3d_hires(datdir, 'QVAPOR', t0, t1, mask=True, drop=True) # kg/kg
    tmpk = var_read_3d_hires(datdir, 'T', t0, t1, mask=True, drop=True) # K
    lw = var_read_3d_hires(datdir, 'RTHRATLW', t0, t1, mask=True, drop=True) # K/s
    lwc = var_read_3d_hires(datdir, 'RTHRATLWC', t0, t1, mask=True, drop=True) # K/s
    sw = var_read_3d_hires(datdir, 'RTHRATSW', t0, t1, mask=True, drop=True) # K/s
    swc = var_read_3d_hires(datdir, 'RTHRATSWC', t0, t1, mask=True, drop=True) # K/s

    print("Computing means")
    w_imean, rho_imean, qv_imean, tmpk_imean, \
    lw_imean, lwc_imean, sw_imean, swc_imean = compute_mean_profiles(nmean, pclass, cwv, w, rho, qv, tmpk, lw, lwc, sw, swc)

    # Save variables for each ens member
    # w_mean[itest,imemb,...]    = w_imean
    # qv_mean[itest,imemb,...]   = qv_imean
    # tmpk_mean[itest,imemb,...] = tmpk_imean
    w_mean[itest,...]    = w_imean
    rho_mean[itest,...]  = rho_imean
    qv_mean[itest,...]   = qv_imean
    tmpk_mean[itest,...] = tmpk_imean
    lw_mean[itest,...]   = lw_imean
    lwc_mean[itest,...]  = lwc_imean
    sw_mean[itest,...]   = sw_imean
    swc_mean[itest,...]  = swc_imean

# ### Save processed data as pickle

# if do_write:
# if comm.rank > 0:

    # comm.Send(np.ascontiguousarray(w_mean, dtype=np.float64), dest=0, tag=comm.rank+20)
    # comm.Send(np.ascontiguousarray(qv_mean, dtype=np.float64), dest=0, tag=comm.rank+30)
    # comm.Send(np.ascontiguousarray(tmpk_mean, dtype=np.float64), dest=0, tag=comm.rank+40)

# if comm.rank == 0:

    # comm.Send(np.ascontiguousarray(w_mean, dtype=np.float64), dest=0, tag=comm.rank+20)
    # comm.Send(np.ascontiguousarray(qv_mean, dtype=np.float64), dest=0, tag=comm.rank+30)
    # comm.Send(np.ascontiguousarray(tmpk_mean, dtype=np.float64), dest=0, tag=comm.rank+40)

pickle_file = main_pickle+memb_all[imemb]+'/mean_profiles_'+str(t1_test)+'hrs.pkl'
with open(pickle_file, 'wb') as file:
    pickle.dump([w_mean, rho_mean, qv_mean, tmpk_mean,
                 lw_mean, lwc_mean, sw_mean, swc_mean], file)
# else:
#     with open(pickle_file, 'rb') as file:
#         w_mean,qv_mean,tmpk_mean = pickle.load(file)

print(comm.rank, 'DONE!!')
