# Friends-of-Friends tools

This folder contains scrip intended to modify or generate Swift FoF group particle
data.

## `remove_fof_groups.py`

The script will remove FoF groups from the particle data if they are below a chosen
size threshold. It is provided for cases when the FoF group minimum size (`FOF:min_group_size`)
was set to be lower than intended, which directly affects the subhalo population that
HBT-HERONS finds due to its reliance on a FOF catalogue.

Note that the old FOF group IDs of each particle are backed up (`FOFGroupIDs -> FOFGroupIDs_old`)
in case something goes wrong during the execution of this script.

If many snapshots require this correction, we recommend using a SLURM array job.
