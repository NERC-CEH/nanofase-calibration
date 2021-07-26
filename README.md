# NanoFASE model calibration

This repository is a template for calibrating the [NanoFASE model](https://github.com/nerc-ceh/nanofase). It is intended for expert use and only as an example of how you may wish to calibrate the model, not as the definitive method of doing so.

The Python file at [`scripts/optimize.py`](./scripts/optimize.py) is an example script of how calibration of the model can be performed, and you may wish to copy this script and modify it for you own needs. Example data is provided in the `data/` and `obs_data/` directories, which is for the Thames catchment between 2009 and 2012.

## Usage

To use the example [optimize.py](./scripts/optimize.py) script with the example data provided in `data/` and `obs_data`, follow these instructions. Firstly, let's clone this repo and make some directories:

```bash
$ git clone https://github.com/nerc-ceh/nanofase-calibration
$ cd nanofase-calibration
$ mkdir config_cache data_cache run_stdout results output
```
Make sure the model has been compiled and that an executable resides somewhere, let's say at `/path/to/model/exe`.

The Python packages required to run the optimize.py script are listed in the Conda [environment.yaml](./environment.yaml) file. Either create a new Conda environment using this:

```bash
$ conda env create -f environment.yaml
$ conda activate nanofase-calibration
```

Or install the required packages however else you wish. Now you can run the optimize.py script, passing in the path to the calibration directory and the model executable:

```bash
(nanofase-calibration) $ cd scripts
(nanofase-calibration) $ ./optimize.py --caldir ../ --exepath /path/to/model/exe
```

When the calibration has finished, the final calibration parameters are stored in `results/optimized_params.npy`, which is a NumPy binary file. To read this:

```bash
(nanofase-calibration) $ python
>>> import numpy as np
>>> params = np.load('results/optimized_params.npy')
```

### Saving intermediate results

To save memory, the optimization script deletes input and output data after each iteration has completed. This means that, if the script crashes, you won't have any information on how far the iteration got. To avoid this, you can specify the `--logall` flag, which results in a `.npz` being written to the `results/` directory for each iteration (of which there could be 1000s), with the run ID hash as the file name. The disadvantage of this is that these files could take up a significant amount of disk space. To inspect these, use NumPy:

```bash
(nanofase-calibration) $ python
>>> import numpy as np
>>> result = np.load('results/<run_id_hash>.npz')
>>> cost = result['cost']
>>> params = result['params']
```

The `cost` variable above is the mean absolute error compute for that particular model run.

Use the optimize.py script's `-h` flag to see all options:

```
usage: optimize.py [-h] [--caldir CALDIR] [--exepath EXEPATH] [--logall] [--no-logall] [--maxworkers MAXWORKERS]
                   [--yearrange YEARRANGE YEARRANGE]

Optimize the NanoFASE model

optional arguments:
  -h, --help            show this help message and exit
  --caldir CALDIR, -c CALDIR
                        path to the calibration directory
  --exepath EXEPATH, -x EXEPATH
                        path to the model exe file
  --logall, -l
  --no-logall
  --maxworkers MAXWORKERS, -mw MAXWORKERS
                        maximum number of processes that can be used
  --yearrange YEARRANGE YEARRANGE, -yr YEARRANGE YEARRANGE
                        year range to run calibration for (inclusive)
```

### Cleaning up

The [`clean`](./clean) bash script is provided to clean up results and the cache from a calibration run.

## How the calibration works

The model calibration is done on 8 parameters, which control suspended sediment concentrations. Each of these parameters is allowed to vary spatially (but not temporally), and so the total number of parameters that are being calibrated is 8 multiplied by the number of grid cells in your geographical scenario (which, as you've probably guessed, could amount to *a lot* of parameters). Suspended sediment observation data at given sampling sites is used to calibrate the model against. The general goal of the calibration is to minimise the *mean absolute error* of the simulated data vs the observation data, and this minimisation is done in parallel by the [optimparallel](https://pypi.org/project/optimparallel/) package. Internally, this package uses the [L-BFGS-B method](https://en.wikipedia.org/wiki/Limited-memory_BFGS), thereby providing a parallel version of `scipy.optimize.minimize(method='L-BFGS-B')`.

Have a look at the example observation data in `obs_data/` to see what format the observation is required in. Note that a [`sites_meta.csv`](./obs_data/sites_meta.csv) file is required to provide metadata (e.g. grid references), whilst the observation data itself is provided in [`sites_obs.csv`](./obs_data/sites_obs.csv).

Initial values for these parameters should be provided in the `data/x0.nc` file.

## Notes and caveats

- Depending on your model scenario (e.g. spatial and temporal extent), the calibration could take *a very long time*. For a desktop computer calibrating the Thames catchment for 2009 to 2012, we could be talking months or years. Two approaches are suggested to alleviate this:
	- The optimize.py script assumes you want to vary each calibration parameter (of which there are 8) for each grid cell (of which, in the Thames, there are 634). Instead, you could assume that the calibration parameters are constant over the whole catchment, or split the catchment into subcatchments and assume the calibration parameters are constant over each of these.
    - Find a high-performance computer with a ridiculous number of cores (>100s). The optimize.py script automatically uses the maximum number of cores that it can (or you can specify the maximum number to use by passing specifying the `--maxworkers` argument). It will still probably take a while.
- Make sure the compiled model executable you are using was compiled with speed optimization in mind. If you are using the example Makefile, then use `make release` or `make fast` to achieve this. This can more than halve model run times.
- The optimize.py script has only been tested using Bash and may need modifications to use on different systems, particularly Windows.