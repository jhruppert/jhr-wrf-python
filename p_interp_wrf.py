# Purpose: Vertically interpolate your WRF datasets!

# Input:
    # input_file: The stitched WRFout file path
    # variable_name: A list of variables that you are interested in calculating
        # U == Zonal wind [m/s]
    # output_dir: The path to a directory where you'd like the new .nc files to be located
    # vertical_levels: An np.array() of pressure level(s) (in hPa) to interpolat
# Output:
    # .nc files for specific variables
# Process:
    # Open the stitched wrfout file
    # Figure out if the user wants 1 level or multiple levels, then loop through the variables
    # Create the new .nc file, copy global attributes over, and edit certain dims
    # Create the home where the variable will live then loop through each timestep
        # and fill it with the interpolated variable. This loop is necessary for 
        # variables that are too big to load into one variable.
# Tip:
    # You'd want to run this function for each domain file you have because input_file currently takes one path.
## EXAMPLE ##
# i.e. if I want to interpolate zonal winds on pressure coordinates on 50hPa , I would run this: 
# parent_dir = '/this/is/where/my/data/lives'
# input_file_d01 = parent_dir + '/raw/d01'  # Path to the raw input netCDF file
# input_file_d02 = parent_dir + '/raw/d02'  # Path to the raw input netCDF file
# output_dir = parent_dir + '/L2/'  # Path to the directory with interpolated files
# variable_name = ['U']             # Declare the variables you want to interpolate
# vertical_levels = np.arange(1000,0,-50)   # Pressure heights you want to interpolate at
# Call the function:
# interp_variable(input_file_d01, variable_name, output_dir, vertical_levels)
# 
# Hrag Najarian, University of Oklahoma
# 1 June 2023

import netCDF4 as nc
import numpy as np
import wrf

def interp_variable(input_file, variable_name, output_dir, vertical_levels):
    # Open the input netCDF file
    dataset = nc.Dataset(input_file, 'r')   # 'r' is just to read the dataset, we do NOT want write privledges

    if vertical_levels.shape == (): levels = 1
    else: levels = len(vertical_levels)

    for i in variable_name:
        if i == 'U':
            # Create new .nc file we can write to and name it appropriately
            if levels == 1:
                output_dataset = nc.Dataset(output_dir + input_file[-3:] + '_interp_' + 'U' + str(vertical_levels), 'w', clobber=True)
            else:
                output_dataset = nc.Dataset(output_dir + input_file[-3:] + '_interp_U', 'w', clobber=True)
            output_dataset.setncatts(dataset.__dict__)
            # Create the dimensions based on global dimensions, with exception to bottom_top
            for dim_name, dim in dataset.dimensions.items():
                if dim_name == 'bottom_top':    output_dataset.createDimension(dim_name, levels)
                else:   output_dataset.createDimension(dim_name, len(dim))
            # wrf.getvar() will destagger, therefore dim 'west_east_stag' should be 'west_east'
            temp = list(dataset.variables['U'].dimensions)
            temp[-1] = 'west_east'
            temp = tuple(temp)
            # Create the variable, set attributes, and start filling the variable into the new nc file
            output_variable = output_dataset.createVariable(i, 'f4', temp)  # 'f4' == float32
            temp_atts = dataset.variables['U'].__dict__
            temp_atts.update({'stagger': ''})
            output_variable.setncatts(temp_atts)
            for t in range(dataset.dimensions['Time'].size):
                variable = wrf.getvar(dataset, 'ua', timeidx=t, meta=False)
                pressure_heights = wrf.getvar(dataset, 'pressure', timeidx=t, meta=False)
                interp_variable = wrf.interplevel(variable,pressure_heights,vertical_levels,meta=False)
                output_variable[t,...] = interp_variable[:]
            # Make sure you close the input and output files at the end
            output_dataset.close()
    dataset.close()
    return