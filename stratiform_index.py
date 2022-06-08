#!/usr/bin/env python
# coding: utf-8

# New simple stratiform index based on rainfall rate and vertical mass flux
# 
# James Ruppert  
# jruppert@ou.edu  
# June 2022


import numpy as np

# Inputs:
#   rain - rainfall rate [mm/hr]
#   vmf  - vertical mass flux integrated over troposphere [kg/m/s]
#   verbose - set to 1 if you want statistics printed, 0 otherwise
# Returns:
#   strat - conv/stratiform index, where 0 = non-raining, 1 = conv, 2 = strat


def stratiform_index(rain,vmf_lower,vmf_upper,verbose):
    
    shape = np.shape(rain)
    strat = np.zeros(shape)
    
    #### Rainfall rate threshold
    # 1. Array initialized as non-raining (0)
    # 2. Stratiform (2): 1 mm/hr > rainfall rate > 0.01 mm/hr
    # 3. Convective (1): rainfall rate ≥ 1 mm/hr
    
    thresh_raining = 1e-1 # mm/hr
    thresh_conv    = 1    # mm/hr
    
    # strat[ ( (rain > thresh_raining) & (rain < thresh_conv) ) ] = 2
    strat[ ( (rain > thresh_raining) ) ] = 1
    # strat[ (rain >= thresh_conv) ] = 1
    
    ncell=shape[0]
    for idim in range(1,rain.ndim):
        ncell*=shape[idim]
    
    if verbose == 1:
        
        print("After rainfall threshold:")

        indices = (strat == 0).nonzero()
        count = np.shape(indices[1])[0]
        print("N non-raining = ",count)
        print("% non-raining = ",count/ncell*100)
        print()

        indices = (strat == 1).nonzero()
        count = np.shape(indices[1])[0]
        print("N conv = ",count)
        print("% conv = ",count/ncell*100)
        print()

        indices = (strat == 2).nonzero()
        count = np.shape(indices[1])[0]
        print("N strat = ",count)
        print("% strat = ",count/ncell*100)
        print()
    
    
    #### VMF threshold
    # 1. 1 changed to 2 (C-->S): VMF < 10**3 kg/m/s
    # 2. 2 changed to 0 (S-->Nr): VMF < 10 kg/m/s
    
    vmf = vmf_lower + vmf_upper
    
    thresh_conv   = 1e3  # kg/m/s
    thresh_nonrain = 10   # kg/m/s

    # strat[ ( (strat == 1) & (vmf < thresh_strat) ) ] = 2
    # strat[ ( (strat == 2) & (vmf < thresh_nonrain) ) ] = 0
    # strat[ (vmf >= thresh_nonrain) ] = 2
    # strat[ ((vmf >= thresh_conv) ) ] = 1

    if verbose == 1:
        
        print("After VMF threshold:")
        
        indices = (strat == 0).nonzero()
        count = np.shape(indices[1])[0]
        print("N non-raining = ",count)
        print("% non-raining = ",count/ncell*100)
        print()

        indices = (strat == 1).nonzero()
        count = np.shape(indices[1])[0]
        print("N conv = ",count)
        print("% conv = ",count/ncell*100)
        print()

        indices = (strat == 2).nonzero()
        count = np.shape(indices[1])[0]
        print("N strat = ",count)
        print("% strat = ",count/ncell*100)
        print()
    
    #### VMF TOPHEAVINESS threshold
    # 1. TH = (VMF_up - VMF_low)/VMF_tot
    # 2. 1 changed to 2 (C-->S): TH > 0.75
    
    th = (vmf_upper - vmf_lower) / (vmf_upper + vmf_lower)
    
    thresh_strat   = .75 # %
    
    strat[ ( (strat == 1) & (th > thresh_strat) ) ] = 2

    if verbose == 1:
        
        print("After VMF TOPHEAVINESS threshold:")
        
        indices = (strat == 0).nonzero()
        count = np.shape(indices[1])[0]
        print("N non-raining = ",count)
        print("% non-raining = ",count/ncell*100)
        print()

        indices = (strat == 1).nonzero()
        count = np.shape(indices[1])[0]
        print("N conv = ",count)
        print("% conv = ",count/ncell*100)
        print()

        indices = (strat == 2).nonzero()
        count = np.shape(indices[1])[0]
        print("N strat = ",count)
        print("% strat = ",count/ncell*100)
        print()
        
    
    return strat

