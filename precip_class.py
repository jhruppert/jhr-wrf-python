# 
# Column-based precipitation classification algorithm for numerical model output.
# 
# Input:
# 
#       Q_VAR: 4D array as f(q, z, x1, x2), where q(5) is the hydrometeor
#               dimension, arranged as ['QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP']
# 
#       Q_INT: 3D array as f(q, x1, x2), vertically integrated Q_VAR.
# 
# Returns:
# 
#       C_TYPE: 2D array as f(x1, x2) with classification results:
#               0: non-precipitating
#           Convective:
#               1: deep convective
#               2: congestus
#               3: shallow
#           Layered:
#               4: stratiform
#               5: anvil (weaker rainfall)
# 
# Emily Luschen - emily.w.luschen-1@ou.edu
# James Ruppert - jruppert@ou.edu
# 5/18/23

import numpy as np

def precip_class(q_var, q_int):

    nx1,nx2 = q_var.shape[2:4]

    # Integrated water variables
    LWP = q_int[0] + q_int[1]               # Liquid water path = cloud + rain
    IWP = q_int[2] + q_int[3] + q_int[4]    # Ice water path = ice + snow + graupel
    TWP = LWP + IWP                         # Total water path [mm]

    # Threshold p]arameters
    twp_thresh = 1e-1
    cr_thresh = 2
    # ice_thresh = 1e-8
    graup_thresh = 1e-4
    rain_thresh_conv = 1e-1
    rain_thresh_strat = 1e-2

    # Initialize output array
    c_type = np.zeros((nx1,nx2))

    for ix1 in range(nx1): # loop through lat
        for ix2 in range(nx2): # loop through lon

            # Considered non-cloud/non-precipitating if below TWP threshold
            if TWP[ix1,ix2] > twp_thresh:

                # Skip points with zero CWP (NANs)
                if LWP[ix1,ix2] != 0:

                    CR = IWP[ix1,ix2]/LWP[ix1,ix2] # cloud ratio

                    # Convective types
                    if CR <= cr_thresh: # bottom-heaviness threshold

                        if q_int[1,ix1,ix2] >= rain_thresh_conv: # rain threshold
                            if q_int[4,ix1,ix2] >= graup_thresh: # graupel threshold
                                c_type[ix1,ix2] = 1 # deep convective
                            else:
                                c_type[ix1,ix2] = 2 # congestus convective
                        else:
                            c_type[ix1,ix2] = 3 # shallow

                    # Stratiform types
                    else:
                        if q_int[1,ix1,ix2] >= rain_thresh_strat: # rain threshold
                            c_type[ix1,ix2] = 4 # stratiform
                        else:
                            c_type[ix1,ix2] = 5 # anvil

    return c_type