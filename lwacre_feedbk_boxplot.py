# ### Notebook to genereate boxplots for ddt(CWV) binned by cloud classification scheme.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 12/15/23

import numpy as np
import matplotlib
from matplotlib import ticker
import matplotlib.pyplot as plt
import sys
from thermo_functions import *
from precip_class import *
import seaborn as sns
import xarray as xr
from memory_usage import *
from read_functions import *

# #### Main settings

storm = 'haiyan'
# storm = 'maria'

# main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrfenkf/"
main = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/"
datdir2 = 'post/d02/'

# Tests to read and compare
# tests = ['crfon','ncrf']
# if storm == 'haiyan':
#     tests = ['ctl','ncrf36h']
# elif storm == 'maria':
#     # tests = ['ctl','ncrf36h']
#     tests = ['ctl','ncrf48h']
# tests = ['ctl']
test_str='ctl'#tests[ktest]

figdir = "/home/jamesrup/figures/tc/ens/boxplot/"

time_neglect=12 # time steps from start to neglect

# Number of sample time steps
nt=200 # will be chopped down to max available
nt=48

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
datdir = main+storm+'/'+memb_all[0]+'/'+test_str+'/'+datdir2
nt_data, nz, nx1, nx2, pres = get_file_dims(datdir)
dp = (pres[1]-pres[0])*1e2 # Pa
nt=np.min([nt,nt_data-time_neglect])
nx1-=80*2
nx2-=80*2


# Main read loops for 3D (dependent) variables

# Arrays to save variables
# ntest=len(tests)
dims = (nmem, nt, nx1, nx2)
strat_all  = np.ma.zeros(dims)
cwv_1hr    = np.ma.zeros(dims)
lwacre_1hr = np.ma.zeros(dims)
qrain_1hr  = np.ma.zeros(dims)
# vmfu_1hr  = np.ma.zeros(dims)
# condh_1hr  = np.ma.zeros(dims)

# This has been tested for corresponding time steps:
#   t0=37,1 are the first divergent time steps in CTL,NCRF
#   t0=25,1 are the first divergent time steps in NCRF,CRFON
# if test_str == 'ctl':
#     if tests[1] == 'ncrf36h':
#         t0=36
#     elif tests[1] == 'ncrf48h':
#         t0=48
# elif test_str == 'ncrf36h':
#     t0=t0_test
# elif test_str == 'ncrf48h':
#     t0=t0_test
# elif test_str == 'crfon':
#     t0=0
t0=time_neglect # neglect the first 12 time steps
t1=t0+nt

# t0+=1 # add one time step since NCRF(t=0) = CTL
# t1 = t0+nt

print('Running test: ',test_str)

# Loop over ensemble members

for imemb in range(nmem):

    print('Running imemb: ',memb_all[imemb])

    datdir = main+storm+'/'+memb_all[imemb]+'/'+test_str+'/'+datdir2
    print(datdir)

    # Stratiform ID
    q_int = read_qcloud(datdir,t0,t1,drop=True) # mm
    strat = precip_class(q_int)

    # CWV
    varname='PW'
    cwv = var_read_2d(datdir,varname,t0,t1,drop=True) # mm
    # ddtq = np.gradient(lwnet, axis=0) # mm/hr

    # Save variables for each ens member
    strat_all[imemb,:,:,:]  = strat
    qrain_1hr[imemb,:,:,:]  = q_int[1] # mm
    cwv_1hr[imemb,:,:,:]    = cwv
    lwacre_1hr[imemb,:,:,:] = read_lwacre(datdir,t0,t1,drop=True) # W/m2
    # vmfu_1hr[imemb,:,:,:]   = read_mse_diag(datdir,'vmfu',2,t0,t1,drop=True) # kg/m/s
    # condh_1hr[imemb,:,:,:]  = read_mse_diag(datdir,'condh',2,t0,t1,drop=True) # mm/day


# Smoothing function
# def time_smooth_var(var, nwindow):
#     data_xr = xr.DataArray(var,
#                             coords={'test':np.arange(ntest), 'memb':memb_nums, 'time':np.arange(nt),
#                             'y':np.arange(nx1), 'x':np.arange(nx2)},
#                             dims=['test','memb','time','y','x'])
#     data_smooth = data_xr.rolling(time=nwindow, center=True).mean()
#     return data_smooth.to_masked_array()

# ### Run binning

def get_kstrat_cells(var_in, strat):
    var_indexed = []
    nstrat=6
    for kstrat in range(nstrat):
        indices = (strat == kstrat).nonzero()
        indexed_var = var_in[indices[0],indices[1],indices[2],indices[3]]
        var_indexed.append(indexed_var)
    return var_indexed

qrain_indexed_1hr    = get_kstrat_cells(qrain_1hr, strat_all)
cwv_indexed_1hr      = get_kstrat_cells(cwv_1hr, strat_all)
lwacre_indexed_1hr   = get_kstrat_cells(lwacre_1hr, strat_all)
# vmfu_indexed_1hr     = get_kstrat_cells(vmfu_1hr, strat_all)
# condh_indexed_1hr    = get_kstrat_cells(condh_1hr, strat_all)
# lwfeedb1_indexed_1hr = get_kstrat_cells(np.absolute(lwacre_1hr)/qrain_1hr, strat_all)
# lwfeedb2_indexed_1hr = get_kstrat_cells(np.absolute(lwacre_1hr)/vmfu_1hr , strat_all)
# lwfeedb3_indexed_1hr = get_kstrat_cells(np.absolute(lwacre_1hr)/condh_1hr, strat_all)
lwfeedb1_indexed_1hr = get_kstrat_cells(lwacre_1hr/qrain_1hr, strat_all)
# lwfeedb2_indexed_1hr = get_kstrat_cells(lwacre_1hr/vmfu_1hr , strat_all)
# lwfeedb3_indexed_1hr = get_kstrat_cells(lwacre_1hr/condh_1hr, strat_all)

# Mask out where values go infinite
nvars = np.size(lwfeedb1_indexed_1hr)
for ivar in range(nvars):
    # lwfeedb1_indexed_1hr[ivar] = np.ma.masked_where(np.bitwise_or((lwacre_indexed_1hr[ivar] == 0), (qrain_indexed_1hr[ivar] == 0)),
    #     lwfeedb1_indexed_1hr[ivar], copy=False)
    lwfeedb1_indexed_1hr[ivar] = np.ma.masked_where((qrain_indexed_1hr[ivar] == 0), lwfeedb1_indexed_1hr[ivar], copy=False)

# # Mask out where values go infinite
# for ivar in range(nvars):
#     lwfeedb2_indexed_1hr[ivar] = np.ma.masked_where(np.bitwise_or((lwacre_indexed_1hr[ivar] == 0), (vmfu_indexed_1hr[ivar] == 0)),
#         lwfeedb2_indexed_1hr[ivar], copy=False)

# # Mask out where values go infinite
# for ivar in range(nvars):
#     lwfeedb2_indexed_1hr[ivar] = np.ma.masked_where(np.bitwise_or((lwacre_indexed_1hr[ivar] == 0), (condh_indexed_1hr[ivar] <= 0)),
#         lwfeedb2_indexed_1hr[ivar], copy=False)

# ---
# ### Plotting routines

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


# Global boxplot settings

def create_boxplot(invar, title_tag, fig_tag, units, yscale='linear'):

    c_name = ['Non-precip', 'Deep\nConvective', 'Congestus', 'Shallow\nConvective', 'Stratiform', 'Anvil']
    cmap = ['white', 'teal', 'plum', 'darkorange', 'gold', 'cornflowerblue']
    # c_name = ['Deep\nConvective', 'Congestus', 'Shallow\nConvective', 'Stratiform', 'Anvil']
    # cmap = ['teal', 'plum', 'darkorange', 'gold', 'cornflowerblue']
    sns.set_palette(cmap)

    fig = plt.figure(figsize=(5.5,4),dpi=300)
    # fig.set_facecolor('white')
    ax = fig.subplots(nrows=1, ncols=1)
    sns.boxplot([invar[0], invar[1], invar[2], invar[3], invar[4], invar[5]],
    # sns.boxplot([invar[1], invar[2], invar[3], invar[4], invar[5]],
                width=0.7, showmeans=True, #log_scale=log_scale,
                meanprops={"marker":"o", "markerfacecolor":"white", 
                "markeredgecolor":"black", "markersize":"6"})

    ax.set_yscale(yscale)
    ax.set_xticklabels(c_name)
    plt.ylabel(units)#, weight='bold')
    plt.title("Class Averaged "+title_tag)#, weight='bold')
    plt.savefig(figdir+fig_tag+test_str+'.png',dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.2)


print("RUNNING CWV")

units = "mm"
title_tag = "CWV"
fig_tag = "cwv"
create_boxplot(cwv_indexed_1hr, title_tag, fig_tag,  units)

print("DONE")
print("RUNNING LWACRE")

units = "W/m$^2$"
title_tag = "LWACRE"
fig_tag = "lwacre"
create_boxplot(lwacre_indexed_1hr, title_tag, fig_tag, units)

print("DONE")
print("RUNNING VMFup")

# units = "kg/m/s"
# title_tag = "VMFup"
# fig_tag = "vmfu"
# create_boxplot(vmfu_indexed_1hr, title_tag, fig_tag, units, yscale='log')

# print("DONE")
# print("RUNNING CONDH")

# title_tag = "CONDH"
# units = "mm/day"
# fig_tag = "condh"
# create_boxplot(condh_indexed_1hr, title_tag, fig_tag, units, yscale='log')


def create_boxplot_noclear(invar, title_tag, fig_tag, units, yscale='linear'):

    # c_name = ['Non-precip', 'Deep\nConvective', 'Congestus', 'Shallow\nConvective', 'Stratiform', 'Anvil']
    # cmap = ['white', 'teal', 'plum', 'darkorange', 'gold', 'cornflowerblue']
    c_name = ['Deep\nConvective', 'Congestus', 'Shallow\nConvective', 'Stratiform', 'Anvil']
    cmap = ['teal', 'plum', 'darkorange', 'gold', 'cornflowerblue']
    sns.set_palette(cmap)

    fig = plt.figure(figsize=(5.5,4),dpi=300)
    # fig.set_facecolor('white')
    ax = fig.subplots(nrows=1, ncols=1)
    # sns.boxplot([invar[0], invar[1], invar[2], invar[3], invar[4], invar[5]],
    sns.boxplot([invar[1], invar[2], invar[3], invar[4], invar[5]],
                width=0.7, showmeans=True, #log_scale=log_scale,
                meanprops={"marker":"o", "markerfacecolor":"white", 
                "markeredgecolor":"black", "markersize":"6"})

    ax.set_yscale(yscale)
    # ax.set_ylim([1e-2,1e14])
    ax.set_xticklabels(c_name)
    plt.ylabel(units)#, weight='bold')
    plt.title("Class Averaged "+title_tag)#, weight='bold')
    plt.savefig(figdir+fig_tag+test_str+'.png',dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.2)


print("DONE")
print("RUNNING QRAIN")

units = "mm"
title_tag = "QRAIN"
fig_tag = "qrain"
create_boxplot_noclear(qrain_indexed_1hr, title_tag, fig_tag, units, yscale="log")

print("DONE")
print("RUNNING LWFB1")

# Local LW Feedback
units = "W/m$^2$ / mm"
title_tag = "LWACRE/QRAIN"
fig_tag = "lwfdb_qrain"
create_boxplot_noclear(lwfeedb1_indexed_1hr, title_tag, fig_tag, units, yscale="log")

print("DONE")
print("RUNNING LWFB2")

# # Local LW Feedback
# units = "W/m$^2$ / kg/m/s"
# title_tag = "LWACRE/VMFup"
# fig_tag = "lwfdb_vmfu"
# create_boxplot_noclear(lwfeedb2_indexed_1hr, title_tag, fig_tag, units, yscale="log")

# print("DONE")
# print("RUNNING QRAIN")

# Local LW Feedback
# units = "W/m$^2$ / mm/day"
# title_tag = "LWACRE/CONDH"
# fig_tag = "lwfdb_condh"
# create_boxplot_noclear(lwfeedb3_indexed_1hr, title_tag, fig_tag, units, yscale="log")
