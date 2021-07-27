#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess
import shutil
import uuid
import pandas as pd
import numpy as np
import f90nml
from netCDF4 import Dataset
from optimparallel import minimize_parallel

# Parse the input arguments
parser = argparse.ArgumentParser(description='Calibrate the NanoFASE model')
parser.add_argument('--caldir', '-c', help='path to the calibration directory', default='./')
parser.add_argument('--exepath', '-x', help='path to the model exe file', default='./nanofase')
parser.add_argument('--logall', '-l', dest='log_all', action='store_true')
parser.add_argument('--no-logall', dest='log_all', action='store_true')
parser.add_argument('--maxworkers', '-mw', help='maximum number of processes that can be used',
                    type=int)
parser.add_argument('--yearrange', '-yr', nargs=2, type=int, help='year range to run calibration for (inclusive)')
parser.add_argument('--test', '-t', choices=['parallel', 'model'])
parser.set_defaults(log_all=False)
parser.set_defaults(yearrange=[2009,2012])
args = parser.parse_args()
cal_dir = args.caldir
exe_path = args.exepath
max_workers = args.maxworkers
year_range = range(args.yearrange[0], args.yearrange[1]+1)

# The parameters we wish to optimize 
param_names = ['resuspension_alpha', 'resuspension_beta', 'sediment_transport_a', 'sediment_transport_c',
               'deposition_alpha', 'deposition_beta', 'bank_erosion_alpha', 'bank_erosion_beta']

# COST FUNCTION
# -------------
# Define the cost function to calibrate using. In this case, using the mean absolute error
def cost(df):
    mae = []
    # Loop through all of the sites in this subcatchment
    for i, site in df_meta.iterrows():
        if site['full-name'] in df_obs['site_name'].str.upper().values:
            # Get the obs date for this site
            df_obs_site = df_obs[df_obs['site_name'].str.upper() == site['full-name']]
            # Slice the output data to get only sim data for this site
            df_sim = df[(df.x == site.x) & (df.y == site.y) & (df.w == site.w)]
            # Merge the sim and obs data
            df_merged = pd.merge(left=df_obs_site, right=df_sim, left_on='date', right_on='datetime')
            # Calculate the RMS for this site
            with np.errstate(divide='ignore'):
                arr_C_spm = np.log(df_merged['C_spm(kg/m3)']*1e3).replace(-np.Inf, 0)
                df_merged['err'] = abs(arr_C_spm - np.log(df_merged['sld_sus(mg/l)']))
            df_merged['err'] = df_merged['err'].replace(-np.Inf, 0)
            mae.append(df_merged['err'].mean())
    # Average the mean absolute errors of each site
    cost = np.array(mae).mean()
    return cost

# Read in the data to calibrate using, including site metadata
df_meta = pd.read_csv(os.path.join(cal_dir, 'obs_data/sites_meta.csv'))
df_obs = pd.read_csv(os.path.join(cal_dir, 'obs_data/sites_obs.csv'), parse_dates=['datetime', 'date'])

# Create new config files for this run
with open(os.path.join(cal_dir, 'batch_config.nml'), 'r') as f:
    bc_template = f90nml.read(f)
with open(os.path.join(cal_dir, 'config.nml'), 'r') as f:
    config_template = f90nml.read(f)

# Get a mask for the catchment shape and count the cells
nc_subcatchment = Dataset(os.path.join(cal_dir, 'data', f'{args.yearrange[0]}_no-emissions.nc'), 'r')
var = nc_subcatchment['flow_dir'][:,:]
catchment_mask = var.mask
catchment_shape = var.shape
n_cells = var.count()


# MODEL FUNCTION
# --------------
# Function to run the NanoFASE model. This is the function passed to the optimisation function
# and takes the parameters to calibrate as a parameter, returning the cost after running the
# model and analysing the results
def nf_model(params, test=None):
    # Unique run ID to keep track of parallel runs
    run_id = uuid.uuid4().hex
    
    # Log this run
    print(f'Running the model for run ID {run_id}')
    with open(os.path.join(cal_dir, 'optimize.log'), 'a') as f:
        f.write(f'Running the model for run ID {run_id}\n')
    
    # Copy the config templates
    bc = bc_template.copy()
    bc['chunks']['input_files'] = []
    config = config_template.copy()
    config['run']['output_hash'] = run_id
    # Write the config file
    config_path = os.path.join(cal_dir, f'config_cache/config_{run_id}.nml')
    with open(config_path, 'w') as f:
        config.write(f)
    
    # Penalise this iteration if any params less than 0
    if np.any(params < 0.0):
        return 1000
    
    # Make a copy of the template NetCDFs to add this iteration's params to
    for year in year_range:
        dst_path = os.path.join(cal_dir, f'data_cache/{year}_no-emissions_{run_id}.nc')
        shutil.copy(os.path.join(cal_dir, f'data/{year}_no-emissions.nc'), dst_path)
        # Update the batch config to point to this file
        bc['chunks']['input_files'].append(dst_path)
        
    # Write the batch config file
    bc_path = os.path.join(cal_dir, f'config_cache/batch_config_{run_id}.nml')
    with open(bc_path, 'w') as f:
        bc.write(f)
    
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
            
    # Only run model if we're not testing
    if test is None:
        # Run the model and save the output to the run_stdout dir
        with open(os.path.join(cal_dir, f'run_stdout/run_{run_id}.out'), 'w') as file:
            subprocess.run([
                exe_path,
                config_path,
                bc_path,
            ], check=True, text=True, stdout=file, stderr=file)

        # # Remove the output file - we only need it if there was an error, and 
        # # an error will have already trigerred an exception
        os.remove(os.path.join(cal_dir, f'run_stdout/run_{run_id}.out'))
        
        # # Evaluate the output and return the cost
        df_out = pd.read_csv(os.path.join(cal_dir, f'output/output_water{run_id}.csv'),
                            parse_dates=['datetime'])
        # Return the cost, calculated from obs and sim data
        cost_ = cost(df_out)

    # If this is a test, return an arbitrary cost without running the model
    else:
        cost_ = 42

    # If we're logging the params and costs of all iterations, then do so 
    if args.log_all:
        with open(os.path.join(cal_dir, f'results/{run_id}.npz'), 'wb') as f:
            np.savez_compressed(f, cost=cost_, params=params0)

    # Remove the NetCDF files for this run
    for year in year_range:
        os.remove(os.path.join(cal_dir, f'data_cache/{year}_no-emissions_{run_id}.nc'))

    # And the config
    os.remove(os.path.join(cal_dir, f'config_cache/config_{run_id}.nml'))
    os.remove(os.path.join(cal_dir, f'config_cache/batch_config_{run_id}.nml'))
            
    # Also remove the output file - we can recreate this once we've got the calibrated params
    if os.path.isfile(os.path.join(cal_dir, f'output/output_water{run_id}.csv')):
        os.remove(os.path.join(cal_dir, f'output/output_water{run_id}.csv'))
        os.remove(os.path.join(cal_dir, f'output/output_sediment{run_id}.csv'))
        os.remove(os.path.join(cal_dir, f'output/output_soil{run_id}.csv'))
        os.remove(os.path.join(cal_dir, f'output/summary{run_id}.md'))
    
    with open(os.path.join(cal_dir, 'optimize.log'), 'a') as f:
        f.write(f'Cost for {run_id}: {cost_}\n')
    return cost_


# OPTIMIZATION
# ------------
# Perform the optimization, using the x0.nc NetCDF file for initial guesses at
# the parameters
data0 = os.path.join(cal_dir, 'data/x0.nc')
params0 = np.empty((0,))
# For each parameter, get it from the NetCDF file, flatten it and then add to the params0 array
for param in param_names:
    nc = Dataset(data0)
    var = nc[param][:,:]
    flat = np.ravel(var[~var.mask])
    params0 = np.concatenate((params0, flat))

# If this isn't a test, do the optimisation
if args.test in [None, 'parallel']:
    # Actually do the optimization
    result = minimize_parallel(fun=nf_model, x0=params0, args=args.test,
                               parallel={'max_workers': max_workers, 'verbose': True})

    # Write the result to the optimize log
    with open(os.path.join(cal_dir, 'optimize.log'), 'a') as f:
        f.write(f'Final cost (MAE): {result.fun}\n')
        f.write(f'Number of evaluations: {result.nfev}\n')
        f.write(f'Number of iterations: {result.nit}\n')
        f.write(f'Successful? {result.success}')

    # Write params to file, so we can re-run the optimized result
    with open(os.path.join(cal_dir, 'results/optimized_params.npy'), 'wb') as f:
        np.save(f, result.x)

# If this is a model test, just run the model once
elif args.test == 'model':
    cost = nf_model(params0)
    print(cost)
