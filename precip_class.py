# 
# Column-based precipitation classification algorithm designed for application on
# numerical model output.
# 
# It has been designed using WRF model output using the Thompson and Eidhammer
# (2014, JAS) microphysics scheme, which has 2 liquid and 3 frozen categories as
# listed and expected below.
# 
# Input:
# 
#       Q_INT: n-D array of vertically integrated hydrometeors as f(q, X), where
#               q(5) is the hydrometeor dimension, arranged as
#               ['QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP'] and X includes the
#               remaining (time and) spatial dimensions.
# Returns:
# 
#       C_TYPE: (n-2)-D array as f(X) with classification results:
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
# 5/19/23

import numpy as np

def precip_class(q_int):

    shape = q_int.shape
    ndims=len(shape)
    shape_out = shape[1:ndims]

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
    c_type = np.zeros(shape_out)

    cr = IWP/LWP

    # Deep convection
    c_type[( ((LWP != 0) & (TWP > twp_thresh)) &
            (cr <= cr_thresh) &
            (q_int[1] >= rain_thresh_conv) &
            (q_int[4] >= graup_thresh) ).nonzero() ] = 1
    # Congestus
    c_type[( ((LWP != 0) & (TWP > twp_thresh)) &
            (cr <= cr_thresh) &
            (q_int[1] >= rain_thresh_conv) &
            (q_int[4] < graup_thresh) ).nonzero() ] = 2
    # Shallow
    c_type[( ((LWP != 0) & (TWP > twp_thresh)) &
            (cr <= cr_thresh) &
            (q_int[1] < rain_thresh_conv) ).nonzero() ] = 3
    # Stratiform
    c_type[( ((LWP != 0) & (TWP > twp_thresh)) &
            (cr > cr_thresh) &
            (q_int[1] >= rain_thresh_strat) ).nonzero() ] = 4
    # Anvil
    c_type[( ((LWP != 0) & (TWP > twp_thresh)) &
            (cr > cr_thresh) &
            (q_int[1] < rain_thresh_strat) ).nonzero() ] = 5

    return c_type