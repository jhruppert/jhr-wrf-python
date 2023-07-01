# 
# Function to write out a new NetCDF file given a list of variables plus
# their metadata. Note that clobber=True in the Dataset command, so it
# will first delete the file if it exists before writing it.
# 
# Input:
#   - file_out: path+name of NetCDF file to write out
#   - var_list: list containing the variables to write out
#   METADATA:
#   - var_names:  list of string variable names
#   - long_names: list of strings containing variable descriptions
#   - units:      list of strings containing dimensions
#   - dim_names:  list of variable dimensions
# 
# See example code at bottom of script for more info, including an important
# stipulation about the variable dimensions.

from netCDF4 import Dataset
import numpy as np

def write_ncfile(file_out, var_list, var_names, long_names, units, dim_names):

    dims_val = var_list[0].shape

    ncfile = Dataset(file_out,mode='w', clobber=True)

    for idim in range(len(dims_val)):
        dim = ncfile.createDimension(dim_names[0][idim], dims_val[idim]) # unlimited axis (can be appended to).

    nvar = len(var_list)
    for ivar in range(nvar):
        writevar = ncfile.createVariable(var_names[ivar], np.single, dim_names[ivar])
        writevar.units = units[ivar]
        writevar.long_name = long_names[ivar]
        writevar[...] = var_list[ivar]

    ncfile.close()
    
    return None

# # Example code that's been tested

# # A stipulation on the variable list:
# #   All dimensions must appear in the first variable of var_list should,
# #   and any subsequent variables can contain any subcombination of them.
# #   (This could probably be made more flexible...)

# file_out='./test.nc'

# rho  = np.zeros([20, 10, 13, 15])
# rain = np.zeros([20,     13, 15])
# pw   = np.zeros([20,     13, 15])

# var_list=[]
# var_list.append(rho)
# var_list.append(rain)
# var_list.append(pw)

# var_names=['rho','rain','pw']
# long_names=['density of moist air','rainfall rate','precipitable water']
# units=['kg/m^3', 'mm/day', 'mm']
# dim_names=[('nt', 'nz', 'ny', 'nx'), ('nt', 'ny', 'nx'), ('nt', 'ny', 'nx')]

# write_ncfile(file_out, var_list, var_names, long_names, units, dim_names)
