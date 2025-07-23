# Friends-of-Friends tools

This folder contains scrip intended to modify or generate Swift FoF group particle
data.

## `remove_fofs.py`

The script will remove FoF groups from the particle data if they are below a chosen
size threshold. It is provided for cases when the FoF group minimum size (`FOF:min_group_size`)
was set to be lower than intended, which directly affects the subhalo population that
HBT-HERONS finds due to its reliance on a FOF catalogue.

Note that the old FOF group IDs of each particle are backed up (`FOFGroupIDs -> FOFGroupIDs_old`)
in case something goes wrong during the execution of this script.

To use the script, assuming it is being run from this directory, activate the
virtual enviroment and execute the script:
```bash
source ../openmpi-5.0.3-hdf5-1.12.3-env/bin/activate
mpirun -np <NUMBER_MPI_RANKS> -- python3 -m mpi4py remove_fofs.py <SNAP_BASE_PATH> <SNAP_BASE_NAME> <SNAP_IS_DISTRIBUTED> <CATALOGUE_BASE_PATH> <CATALOGUE_IS_DISTRIBUTED> <SNAP_NR> <SIZE_THRESHOLD>
```
If many snapshots require this correction, we recommend using a SLURM array job.