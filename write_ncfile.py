# 
# Function to write out a new NetCDF file given a list of variables plus
# their metadata. Note that clobber=True in the Dataset command, so it
# will first delete the file if it exists before writing it.
# 
# Input:
#   - file_out: path+name of NetCDF file to write out
#   - var_list: list containing the variables to write out
#   METADATA:
#   - var_names:    list of string variable names
#   - descriptions: list of strings containing variable descriptions
#   - units:        list of strings containing units for each variable
#   - dim_names:    list of variable dimensions, i.e., one list or tuple
#                       of dimension strings per variable
# 
# See example code below function for more info, including an important
# stipulation about the variable dimensions.

from netCDF4 import Dataset
import numpy as np

def write_ncfile(file_out, var_list, var_names, descriptions, units, dim_names):

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    dims_val = var_list[0].shape

    ncfile = Dataset(file_out,mode='w', clobber=True)

    for idim in range(len(dims_val)):
        dim = ncfile.createDimension(dim_names[0][idim], dims_val[idim]) # unlimited axis (can be appended to).

    nvar = len(var_list)
    for ivar in range(nvar):
        writevar = ncfile.createVariable(var_names[ivar], np.single, dim_names[ivar])
        writevar.units = units[ivar]
        writevar.description = descriptions[ivar]
        writevar[...] = var_list[ivar]

    ncfile.close()
    
    return None

# Example code that's been tested

# A stipulation on the variable list:
#   All dimensions must appear in the first variable of var_list should,
#   and any subsequent variables can contain any subcombination of them.
#   (This could probably be made more flexible...)

####################################
## Uncomment below #################
####################################

# file_out='./test.nc'

# rho  = np.zeros([20, 10, 13, 15])
# rain = np.zeros([20,     13, 15])
# pw   = np.zeros([20,     13, 15])

# var_list=[]
# var_list.append(rho)
# var_list.append(rain)
# var_list.append(pw)

# var_names=['rho','rain','pw']
# descriptions=['density of moist air','rainfall rate','precipitable water']
# units=['kg/m^3', 'mm/day', 'mm']
# dim_names=[('nt', 'nz', 'ny', 'nx'), ('nt', 'ny', 'nx'), ('nt', 'ny', 'nx')]

# write_ncfile(file_out, var_list, var_names, descriptions, units, dim_names)