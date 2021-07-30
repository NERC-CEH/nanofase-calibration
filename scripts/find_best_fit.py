#!/usr/bin/env python3
import os
import sys
import argparse
import re
import shutil
import numpy as np
from netCDF4 import Dataset

"""This script is useful for finding the best find from the `optimize.log` file,
the getting the parameters for this fit from logged results"""

# Parse the input arguments
parser = argparse.ArgumentParser(description='Find the best parameters so find')
parser.add_argument('--caldir', '-c', help='path to the calibration directory', default='./')
parser.add_argument('--yearrange', '-yr', nargs=2, type=int, help='year range to run calibration for (inclusive)')
parser.set_defaults(yearrange=[2009,2012])
args = parser.parse_args()
cal_dir = args.caldir
year_range = range(args.yearrange[0], args.yearrange[1]+1)

# Open the optimize.log file and find the best fit by plucking the cost and run ID
# from each line that begins with 'C' (Cost for...)
with open(os.path.join(cal_dir, 'optimize.log')) as f:
    costs = []
    ids = []
    for line in f:
        if line[0] == 'C':
            split = re.split('Cost for |\: ', line)
            ids.append(split[1])
            costs.append(float(split[2]))
# Print the minimum cost and the corresponding run ID
costs = np.array(costs)
run_id = ids[costs.argmin()]
print(f'Minimum cost: {costs.min()}')
print(f'For run ID: {run_id}')

# Now get the parameters that produced that cost
params_f = np.load(os.path.join(cal_dir, 'results', f'{run_id}.npz'))
params = params_f['params']

# Finally, we can recreate the NetCDF files used for this run
param_names = ['resuspension_alpha', 'resuspension_beta', 'sediment_transport_a', 'sediment_transport_c',
               'deposition_alpha', 'deposition_beta', 'bank_erosion_alpha', 'bank_erosion_beta']
# Get the template for the 2D array
nc_subcatchment = Dataset(os.path.join(cal_dir, 'data', f'{args.yearrange[0]}_no-emissions.nc'), 'r')
var = nc_subcatchment['flow_dir'][:,:]
catchment_mask = var.mask
catchment_shape = var.shape
n_cells = var.count()
# Make a copy of the template NetCDFs to add this iteration's params to
for year in year_range:
    dst_path = os.path.join(cal_dir, f'data_cache/{year}_no-emissions_{run_id}.nc')
    shutil.copy(os.path.join(cal_dir, f'data/{year}_no-emissions.nc'), dst_path)
# Pull out the 1D arrays for each parameter from the params variable, then
# reshape to the correct grid shape and mask and add to NetCDF file
for i, param in enumerate(param_names):
    param_1d = params[n_cells*i:n_cells*i+n_cells]
    param_2d = np.ma.masked_array(np.empty(catchment_shape), mask=catchment_mask)
    # Reshape into 2D arrays, taking into account the mask
    j = 0
    for i, _ in np.ndenumerate(param_2d):
        if ~catchment_mask[i]:
            param_2d[i] = param_1d[j]
            j = j + 1
    # Now add this variable to the NetCDF file, placing a copy in the cache
    for year in year_range:
        # Then create the new variables
        nc = Dataset(os.path.join(cal_dir, f'data_cache/{year}_no-emissions_{run_id}.nc'), 'r+')
        var = nc.createVariable(param, datatype=float, dimensions=('y','x'))
        var[:] = param_2d
        nc.close()
