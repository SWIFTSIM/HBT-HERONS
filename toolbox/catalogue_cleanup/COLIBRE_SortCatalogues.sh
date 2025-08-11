#!/bin/bash -l
#
# Post-process COLIBRE HBT-HERONS catalogues to join them into a single HDF5 file,
# with the subhaloes sorted in ascending TrackId.
#
# After editing in INDIR & OUTDIR variables below, the script can be run as follows:
#
# mkdir logs
# sbatch -J L0025N0188/Thermal --array=0-127 SortCatalogues.sh
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -o ./logs/SortOutput_%a.%A.out
#SBATCH -p cosma8
#SBATCH -A dp004
#SBATCH -t 01:00:00

set -e

# We assume we are in COSMA
module purge
module load python/3.12.4 gnu_comp/14.1.0 openmpi/5.0.3 parallel_hdf5/1.12.3
source ./swiftsim/openmpi-5.0.3-hdf5-1.12.3-env/bin/activate

# Base input directory where the original HBT catalogues are stored.
INDIR="/cosma8/data/dp004/jlvc76/COLIBRE/ScienceRuns/"

# Base output directory where the joint catalogues will be saved
OUTDIR="/snap8/scratch/dp004/dc-mcgi1/COLIBRE/sort_hbt/"

# Which simulation to do
sim="${SLURM_JOB_NAME}"

# Snapshot index to do
snap_nr=${SLURM_ARRAY_TASK_ID}

# Create ordered catalogue
mpirun -- python -u ./SortCatalogues.py \
 --with-particles \
 --with-potential-energy \
 "${INDIR}/${sim}/HBT-HERONS" \
 "${snap_nr}" \
 "${OUTDIR}/${sim}/HBT-HERONS"

# Set file to read only
snapnum=$(printf "%03d" "$snap_nr")
chmod a=r "${OUTDIR}/${sim}/HBT-HERONS/OrderedSubSnap_${snapnum}.hdf5"

echo "Job complete!"
