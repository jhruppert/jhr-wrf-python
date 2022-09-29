#!/usr/bin/env python
# coding: utf-8
  
# Smooth an array of the form array(time, x, y) or array(time, x, y)
#   using a simple boxcar running mean by leveraging numpy.cumsum
# 
# James Ruppert  
# jruppert@ou.edu  
# June 2022

import numpy as np

#   var_in - assumed to be either of the form:
# 
#                x,y:   array[dim1, dim2]
#                or
#                t,x,y: array[dim1, dim2, dim3]
# 
#   N      - width of boxcar running mean

def smooth_runn_mn(var_in, N):

    # Check for odd
    if N % 2 == 0: N+=1

    shape=np.shape(var_in)

    # If two-dimensional, run one step of smoothing and return
    if len(shape) == 2:
        var_smooth=do_smoothing(ivar_in,N)
        return var_smooth

    # If three-dimensional, run time loop...

    nt=shape[0]
    var_smooth=np.zeros(shape)

    for it in range(nt):
        var_smooth[it,:,:]=do_smoothing(var_in[it,:,:],N)

    return var_smooth



def do_smoothing(ivar_in,N):
    
    shape=np.shape(ivar_in)

    fillval=9e20 # For masking

    ivar_out=np.zeros(shape)
    ivar_out[:]=fillval # Fill for masking

    ind_valid_x1=np.arange((N-1)/2-1, shape[0]-1-(N-1)/2).astype(int)
    ind_valid_x2=np.arange((N-1)/2-1, shape[1]-1-(N-1)/2).astype(int)

    # First dimension
    cumsum = np.cumsum(np.insert(ivar_in, 0, 0, axis=0), axis=0)
    x_sm=(cumsum[N:,:] - cumsum[:-N,:]) / float(N)
    ivar_out[ind_valid_x1,:]=x_sm

    # Second dimension
    cumsum = np.cumsum(np.insert(ivar_in, 0, 0, axis=1), axis=1)
    x_sm=(cumsum[:,N:] - cumsum[:,:-N]) / float(N)
    ivar_out[:,ind_valid_x2]=x_sm

    # Mask out edges
    ivar_out[:np.min(ind_valid_x1)-1,:]=fillval
    ivar_out[-np.min(ind_valid_x1)-1:,:]=fillval
    ivar_out[:,:np.min(ind_valid_x2)-1]=fillval
    ivar_out[:,-np.min(ind_valid_x2)-1:]=fillval

    ivar_out=np.ma.masked_equal(ivar_out,fillval,copy=False)

    return ivar_out


# for ix in range(shape[0]):
#     cumsum = np.cumsum(np.insert(var_in[ix,:], 0, 0)) 
#     x_sm2=(cumsum[N:] - cumsum[:-N]) / float(N)
#     print(np.shape(x_sm2))
#     sys.exit()
#     var_out[ix,ind_valid]=x_sm
