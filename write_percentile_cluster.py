# Notebook to diagnose cluster information based on percentiles of e.g. rainfall, CWV, ...
# 
# James Ruppert  
# jruppert@ou.edu  
# 10/17/24

# ### Main settings

import numpy as np
from read_functions import *
from scipy.ndimage import label, sum as ndi_sum
import pickle
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()

depend_var_tag = ['rainrate', 'PW']
nvar = len(depend_var_tag)

storm = 'haiyan'
# storm = 'maria'

main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"+storm+'/'
datdir2 = 'post/d02/'
main_pickle = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/jruppert/tc_postproc/"+storm+'/'
# main_pickle = "/Users/jamesruppert/code/tc_postproc/"+storm+'/'

# Tests to read and compare
if storm == 'haiyan':
    tests = ['ctl','STRATANVIL_ON','ncrf36h','STRATANVIL_OFF','STRAT_OFF']
    tests_str = ['CTL','CONVOFF','NCRF','STRATANVOFF','STRATOFF']
    # tests = ['ctl','ncrf36h']
    # tests_str = ['CTL','NCRF']
    # tests = ['crfon','ncrf']
elif storm == 'maria':
    # tests = ['ctl','ncrf36h']
    tests = ['ctl','ncrf48h']
    tests_str = ['CTL','NCRF']
ntest=len(tests_str)

# nt=10
nt=200
time_neglect=12 # number of hours to neglect at the beginning of CTL
t1_test=48 # max of 48

do_tests=False
if do_tests:
    nt=t1_test
else:
    tests = ['ctl']
    tests_str = ['CTL']
    ntest=1

# Members
# nmem = 2
nmem = 10 # number of ensemble members (1-5 have NCRF)

# ### Main processing functions

# #### Percentile settings

# List containing percentile values
def get_percentiles():
    # Percentile array
    # nperc = 30
    # percentiles = np.logspace(-3, 2, num=nperc)
    # percentiles = np.array([1,2,3,4,5,6,7,8,9,
    #           10,20,30,40,50,60,70,80,90,
    #           100,200,300,400,500,600,700,800,900,
    #           1000,2000,3000,4000,5000,6000,7000,8000,9000,
    #           10000,20000,30000,40000,50000,60000,70000,80000,90000,])*1e-3
    percentiles = np.arange(2,100,2)
    return percentiles

# #### Variable read functions

# Ensemble member info
memb0=1 # Starting member to read
memb_nums=np.arange(memb0,nmem+memb0,1)
memb_nums_str=memb_nums.astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

def get_dims():
    test_str='ctl'
    datdir = main+'/'+memb_all[0]+'/'+test_str+'/'+datdir2
    # datdir = "/Users/jamesruppert/code/wrfout/"+storm+'/'+memb_all[0]+'/'+test_str+'/'
    nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
    nt_new=np.min([nt,nt_data-time_neglect])
    # dp = (pres[1]-pres[0])*1e2 # Pa
    nx1-=80*2
    nx2-=80*2
    return nt_new, nz, nx1, nx2

def get_read_time_bounds(test_str):
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
    return t0, t1

mask=True; drop=True

# Read independent variables
def get_dependent_variable(var_tag):

    nt, nz, nx1, nx2 = get_dims()

    # Arrays to save variables
    dims = (ntest, nmem, nt, nx1, nx2)
    ivar_depend = np.zeros(dims)

    for itest in range(ntest):
        test_str = tests[itest]
        t0, t1 = get_read_time_bounds(test_str)
        for imemb in range(nmem):
            datdir = main+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
            ivar_depend[itest, imemb,...] = var_read_2d(datdir,var_tag,t0,t1,mask=mask,drop=drop)
    return ivar_depend

def read_variables_imemb(imemb, test_str):

    datdir = datdir = main+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    t0, t1 = get_read_time_bounds(test_str)

    # Stratiform ID
    # q_int = read_qcloud(datdir,t0,t1,mask=mask,drop=drop) # mm
    # pclass = precip_class(q_int)

    # LWACRE
    lwacre = read_lwacre(datdir,t0,t1,mask=mask,drop=drop) # W/m2

    var_str_all = ['PW','rainrate','HFX','LH']
    cwv  = var_read_2d(datdir,var_str_all[0],t0,t1,mask=mask,drop=drop) # mm
    rain = var_read_2d(datdir,var_str_all[1],t0,t1,mask=mask,drop=drop) # mm/d
    sh   = var_read_2d(datdir,var_str_all[2],t0,t1,mask=mask,drop=drop) # W/m2
    lh   = var_read_2d(datdir,var_str_all[3],t0,t1,mask=mask,drop=drop) # W/m2

    var_str_all = ['vmf','condh']
    vmf   = read_mse_diag(datdir,var_str_all[0],t0,t1,mask=mask,drop=drop) # kg/m/s
    condh = read_mse_diag(datdir,var_str_all[1],t0,t1,mask=mask,drop=drop) # kg*K/(m2*s)
    condh *= 1004. # --> W/m^2
    # netCDF metadata is incorrect - it was never converted from native units: kg*K/(m2*s)

    return cwv, rain, lwacre, sh, lh, vmf, condh

# #### Main read/processing loop

for ivar in range(nvar):

    print()
    print("Running ivar: ",depend_var_tag[ivar])

    # Read dependent variables
    depend_var = get_dependent_variable(depend_var_tag[ivar]) # (ntest, nmem, nt, nx1, nx2)

    # Get percentiles
    percentiles = get_percentiles()
    depend_var_perc = np.nanpercentile(depend_var.flatten(), percentiles)
    # ind = np.nonzero(depend_var_perc)
    # depend_var_perc = depend_var_perc[ind]
    # percentiles = percentiles[ind]
    nperc = percentiles.size

    nt, nz, nx1, nx2 = get_dims()

    for itest in range(ntest):

        test_str = tests[itest]
        print()
        print('Running test: ',test_str)

        # Loop over ensemble members
        # for imemb in range(nmem):
        # for imemb in range(1):
        imemb = comm.rank

        print('Running imemb: ',memb_all[imemb])

        # Arrays to save variables
        # dims = (nvar, ntest, nmem, nt, nperc)
        dims = (nt, nperc)

        number = np.full(dims, np.nan)
        mean_size = np.full(dims, np.nan)

        rain_perc   = np.zeros(dims)
        cwv_perc    = np.zeros(dims)
        lwacre_perc = np.zeros(dims)
        sh_perc     = np.zeros(dims)
        lh_perc     = np.zeros(dims)
        vmf_perc    = np.zeros(dims)
        condh_perc  = np.zeros(dims)

        cwv, rain, lwacre, sh, lh, vmf, condh = read_variables_imemb(imemb,  test_str)
        condh *= 1004. # --> W/m2 # netCDF metadata is incorrect - Condensation heating was never converted from native units: kg*K/(m2*s)

        for it in range(nt):
            for iperc in range(nperc):
                labeled_matrix, inum_features = label(depend_var[itest,imemb,it,...] >= depend_var_perc[iperc])
                cluster_sizes = ndi_sum(depend_var[itest,imemb,it,...], labeled_matrix, index=np.arange(1, inum_features + 1))
                number[it, iperc] = inum_features

                # Average variables over clusters
                ind = np.where(labeled_matrix != 0)
                if ind[0].size != 0:
                    mean_size[it, iperc]  = np.mean(cluster_sizes)
                    rain_perc[it,iperc]   = np.mean(rain[it,ind[0],ind[1]])
                    cwv_perc[it,iperc]    = np.mean(cwv[it,ind[0],ind[1]])
                    lwacre_perc[it,iperc] = np.mean(lwacre[it,ind[0],ind[1]])
                    sh_perc[it,iperc]     = np.mean(sh[it,ind[0],ind[1]])
                    lh_perc[it,iperc]     = np.mean(lh[it,ind[0],ind[1]])
                    vmf_perc[it,iperc]    = np.mean(vmf[it,ind[0],ind[1]])
                    condh_perc[it,iperc]  = np.mean(rain[it,ind[0],ind[1]])

        # Write out to pickle file
        pickle_file = main_pickle+memb_all[imemb]+'/percentile_cluster_'+depend_var_tag[ivar]+'_'+test_str+'_'+str(nt)+'hrs.pkl'
        with open(pickle_file, 'wb') as file:
            pickle.dump([depend_var_perc,number,mean_size,
                            rain_perc,cwv_perc,lwacre_perc,sh_perc,lh_perc,vmf_perc,condh_perc],
                            file)

print("Done!!")
