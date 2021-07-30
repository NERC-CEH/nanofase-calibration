#!/usr/bin/env python3
import os
import argparse
import re
import numpy as np

"""This script is useful for finding the best find from the `optimize.log` file,
the getting the parameters for this fit from logged results"""

# Parse the input arguments
parser = argparse.ArgumentParser(description='Find the best parameters so find')
parser.add_argument('--caldir', '-c', help='path to the calibration directory', default='./')
args = parser.parse_args()
cal_dir = args.caldir

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
id = ids[costs.argmin()]
print(f'Minimum cost: {costs.min()}')
print(f'For run ID: {id}')

# Now get the parameters that produced that cost
params_f = np.load(os.path.join(cal_dir, 'results', f'{id}.npz'))
params = print(params_f['params'])

# Now you can split that params array to get each individual parameter,
# and then potentially create a NetCDF file from it