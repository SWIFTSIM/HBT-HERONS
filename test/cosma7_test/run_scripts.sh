#!/bin/bash

set -e

if [ "$#" -ne 2 ]
then
  echo "Usage: ./run_scripts.sh <PATH_TO_HBT_BASE_FOLDER> <SNAPSHOT_INDEX>"
  exit 1
fi

BASE_FOLDER=${1}
SNAPSHOT_INDEX=${2}

# Check if the SNAPSHOT_INDEX is within the allowed bounds.
MAX_SNAPSHOT_INDEX=$(grep 'MaxSnapshotIndex' $BASE_FOLDER/Parameters.log | awk '{print $2}')
if [ $SNAPSHOT_INDEX -gt $MAX_SNAPSHOT_INDEX ];
then
  echo "You are trying to analyse SNAPSHOT_INDEX = $SNAPSHOT_INDEX but the maximum value is $MAX_SNAPSHOT_INDEX".
  exit 1
fi

# This is unnecessary for this test, as it is always the same simulation. However, I am making
# it a more general in case we change the simulation in the future.
SNAPSHOT_FOLDER=$(grep 'HaloPath' $BASE_FOLDER/Parameters.log | awk '{print $2}')
SNAPSHOT_BASENAME=$(grep 'SnapshotFileBase' $BASE_FOLDER/Parameters.log | awk '{print $2}')

# We need to find the correspondence between snapshot number and snapshot index. These are not
# the same if we did a subset of all snapshots, so we check for it.
if grep -q 'SnapshotIdList' $BASE_FOLDER/Parameters.log
then
  # We get the SnapshotIdList values and check if there is more than one space-separated string.
  # If there is only one, then no submsaple has been specified.
  SNAPSHOT_SUBSAMPLE_NR=($(grep 'SnapshotIdList' $BASE_FOLDER/Parameters.log | awk '{ print }'))

  if [ $(echo ${#SNAPSHOT_SUBSAMPLE_NR[@]}) -gt 1 ];
  then
    SNAPSHOT_NR=${SNAPSHOT_SUBSAMPLE_NR[((${SNAPSHOT_INDEX}+1))]}
  else
    SNAPSHOT_NR=${SNAPSHOT_INDEX}
  fi
else
  SNAPSHOT_NR=${SNAPSHOT_INDEX}
fi

# Load the modules. We assume we have created the relevant virtual enviroment in COSMA
module purge
module load python/3.12.4 gnu_comp/14.1.0 openmpi/5.0.3 parallel_hdf5/1.12.3
source ../../toolbox/swiftsim/openmpi-5.0.3-hdf5-1.12.3-env/bin/activate

# Location of the snapshot
SNAPSHOT_PATH="${SNAPSHOT_FOLDER}/${SNAPSHOT_BASENAME}_{snap_nr:04d}.hdf5" # SWIFT snapshots

# Tests whether particles of subgroups are confined to their host group (or are not part of any FOF)
mpirun -np 16 python3 -m mpi4py ./scripts/check_interhost_subhalos.py "${BASE_FOLDER}" $SNAPSHOT_INDEX $SNAPSHOT_NR --snapshot-file="${SNAPSHOT_PATH}";
echo

# Tests whether there are any duplicated particles in the HBT catalogues.
mpirun -np 16 python3 -m mpi4py ./scripts/check_presence_duplicate_particles.py "${BASE_FOLDER}" $SNAPSHOT_INDEX;
echo

# Tests whether the ID of the tracer used for orphaned subgroups remains unchanged.
mpirun -np 16 python3 -m mpi4py ./scripts/check_consistency_particles_orphan_subgroups.py "${BASE_FOLDER}" $SNAPSHOT_INDEX;
echo

# Tests correctness of assigned host groups for orphaned subgroups.
mpirun -np 16 python3 -m mpi4py ./scripts/check_tracing_particles_orphan_subgroups.py "${BASE_FOLDER}" $SNAPSHOT_INDEX $SNAPSHOT_NR --snapshot-file="${snap_files}";
echo

# Tests correctness of assigned host groups for resolved subgroups.
mpirun -np 16 python3 -m mpi4py ./scripts/check_tracing_particles_resolved_subgroups.py "${BASE_FOLDER}" $SNAPSHOT_INDEX $SNAPSHOT_NR --snapshot-file="${snap_files}";
