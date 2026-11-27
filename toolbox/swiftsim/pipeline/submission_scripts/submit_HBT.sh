#!/bin/bash
#SBATCH --ntasks 16
#SBATCH --cpus-per-task=8
#SBATCH -p cosma8
#SBATCH -A dp004
#SBATCH -t 02:00:00

set -e

module purge
module load gnu_comp/14.1.0 hdf5/1.14.4 openmpi/5.0.3

echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

# Path to HBT executable
HBT_executable=CURRENT_PWD/../../../build/HBT

echo "Using configuration present in ${1}" 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --bind-to none "${HBT_executable}" $1 $2 $3

echo "Job complete!"
