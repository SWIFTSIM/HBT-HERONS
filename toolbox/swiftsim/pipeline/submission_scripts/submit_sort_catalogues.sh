#!/bin/bash -l
#
# Post-process the HBT-HERONS catalogues of a single output to join them into a single
# HDF5 file, with the subhaloes sorted in ascending TrackId. Particle IDs and potential
# energies are not copied across. A Subhalos/SOAPIndex dataset is added, giving the index
# of each subhalo in the corresponding SOAP catalogue. This requires reading the SWIFT
# snapshot, whose location is taken from the config.txt file in the HBT-HERONS output
# folder. The sorted catalogues are saved in a sorted_catalogues directory within the
# HBT-HERONS output folder.
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

# Locate the SWIFT snapshot for this output, needed to compute SOAPIndex
SNAPSHOT_FOLDER=$(grep 'HaloPath' "${1}/config.txt" | awk '{print $2}')
SNAPSHOT_BASENAME=$(grep 'SnapshotFileBase' "${1}/config.txt" | awk '{print $2}')
SNAPSHOT_FILE=$(find "${SNAPSHOT_FOLDER}/" -maxdepth 2 \
  -name "$(printf '%s_%04d.hdf5' "$SNAPSHOT_BASENAME" "$SLURM_ARRAY_TASK_ID")" | head -n 1)
if [ -z "${SNAPSHOT_FILE}" ]; then
  echo "Could not find SWIFT snapshot ${SLURM_ARRAY_TASK_ID} in ${SNAPSHOT_FOLDER}. Exiting."
  exit 1
fi

# Run the code
mpirun -- python3 -u -m mpi4py \
  CURRENT_PWD/../../catalogue_cleanup/SortCatalogues.py \
  --swift-cell-file="${SNAPSHOT_FILE}" \
  "${1}" \
  ${SLURM_ARRAY_TASK_ID} \
  "${1}/sorted_catalogues"

echo "Job complete!"
