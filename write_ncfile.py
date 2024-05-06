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
#   - XX dim_names:    [NO LONGER USING] list of variable dimensions, i.e., one list or tuple
#                       of dimension strings per variable
#   - dims_set:     list of paired lists as:
#                       dims_set = dims_set(n-variable)
#                       dims_set[0] = dims_set[0](2, n-dimension)
#                       dims_set[0][0] = tuple of dimension string names
#                       dims_set[0][1] = tuple of dimension values
# 
# See example code below function for more info.
# 
# Author: James Ruppert
# Date: 07/2023

from netCDF4 import Dataset
import numpy as np

def write_ncfile(file_out, var_list, var_names, descriptions, units, dims_set): #, dim_names

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    # dims_val = var_list[0].shape

    ncfile = Dataset(file_out,mode='w', clobber=True)

    # Add dimensions to file
    # Will iterate over entire variable dimension list but will only attempt to add each once
    dim_added = [] # List to track dimensions that have been added already
    for idimset in range(len(dims_set)):
    #     dim = ncfile.createDimension(dim_names[0][idim], dims_val[idim]) # unlimited axis (can be appended to).
        for idim in range(len(dims_set[idimset][0])):
            if dims_set[idimset][0][idim] in dim_added:
                continue
            dim = ncfile.createDimension(dims_set[idimset][0][idim], dims_set[idimset][1][idim]) # unlimited axis (can be appended to).
            dim_added.append(dims_set[idimset][0][idim])

    for ivar in range(len(var_list)):
        print("  Writing var: ",var_names[ivar])
        writevar = ncfile.createVariable(var_names[ivar], np.single, dims_set[ivar][0]) #dim_names[ivar])
        writevar.units = units[ivar]
        writevar.description = descriptions[ivar]
        writevar[...] = var_list[ivar]

    ncfile.close()

    print("Done writing!")
    
    return None

# Example code that's been tested

####################################
## Uncomment below #################
####################################

# file_out='./test.nc'

# nt=20; nz=10; nx1=13; nx2=15
# dims3d=(nt, nz, nx1, nx2)
# dims2d=(nt, nx1, nx2)
# rho  = np.zeros(dims3d)
# rain = np.zeros(dims2d)
# pw   = np.zeros(dims2d)

# var_list=[]
# var_list.append(rho)
# var_list.append(rain)
# var_list.append(pw)

# var_names=['rho','rain','pw']
# descriptions=['density of moist air','rainfall rate','precipitable water']
# units=['kg/m^3', 'mm/day', 'mm']
# dims3d_names=('nt', 'nz', 'nx1', 'nx2')
# dims2d_names=('nt', 'nx1', 'nx2')
# dims_set = [
#     [dims3d,dims3d_names],
#     [dims2d,dims2d_names],
#     [dims2d,dims2d_names],
# ]


# write_ncfile(file_out, var_list, var_names, descriptions, units, dims_set)
