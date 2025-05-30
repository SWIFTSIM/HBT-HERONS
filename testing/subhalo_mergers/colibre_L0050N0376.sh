#!/bin/bash -l
#
#SBATCH --tasks=64
#SBATCH --cpus-per-task=1
#SBATCH -J colibre
#SBATCH -o ./logs/%x_%A.%a.out
#SBATCH -p cosma8-shm2
#SBATCH -A dp004
#SBATCH -t 2:00:00
#

module purge
module load gnu_comp/14.1.0 openmpi
activate failed_mergers

snap_nr=`printf '%03d' ${SLURM_ARRAY_TASK_ID}`
branch="master"

# Compute core position and extent in phase space
mpirun -- python3 -m mpi4py merging_condition.py \
       "/cosma8/data/dp004/colibre/Runs/L0050N0376/Thermal_non_equilibrium/snapshots/colibre_{snap_nr:04d}/colibre_{snap_nr:04d}.{file_nr}.hdf5" \
       "/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/${branch}/membership/membership_{snap_nr:04d}/membership_{snap_nr:04d}.{file_nr}.hdf5" \
       "/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/${branch}/output/" \
       "/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/${branch}/failed_mergers/" \
       ${snap_nr}

# Identify things which should have merged but did not
python3 ./find_unmerged.py \
        "/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/${branch}/failed_mergers/merging_info_${snap_nr}.hdf5" \
        "/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/${branch}/failed_mergers/failed_mergers_${snap_nr}.hdf5"

# Make plots
python3 ./plot_failed_mergers.py ${snap_nr} ${branch}
