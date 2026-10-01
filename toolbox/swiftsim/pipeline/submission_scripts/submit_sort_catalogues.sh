#!/bin/bash -l
#
# Post-process the HBT-HERONS catalogues of a single output to join them into a single
# HDF5 file, with the subhaloes sorted in ascending TrackId. Particle IDs and potential
# energies are not copied across. The sorted catalogues are saved in a sorted_catalogues
# directory within the HBT-HERONS output folder.
#
# Specify the path to the HBT-HERONS output folder, and which snapshot index to do. The
# latter is specified by the array job index. For example, to sort 128 HBT-HERONS outputs:
#
# sbatch --array=0-127 submit_sort_catalogues.sh <HBT_FOLDER>
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -p cosma8
#SBATCH -A dp004
#SBATCH -t 00:10:00

set -e

# This is assuming we run in COSMA
module purge
module load python/3.12.4 gnu_comp/14.1.0 openmpi/5.0.3 parallel_hdf5/1.12.3
source CURRENT_PWD/../openmpi-5.0.3-hdf5-1.12.3-env/bin/activate

# Path to the HBT-HERONS output folder
echo "Sorting catalogues present in ${1}"

# Run the code
mpirun -- python3 -u -m mpi4py \
  CURRENT_PWD/../../catalogue_cleanup/SortCatalogues.py \
  "${1}" \
  ${SLURM_ARRAY_TASK_ID} \
  "${1}/sorted_catalogues"

echo "Job complete!"
