#!/bin/bash

. ../helper_functions.sh
set -e

# This is where the outputs should be placed
BASE_PATH=$1
HBT_FOLDER=$2

# Parse the SWIFT config file
if [ ! -f $BASE_PATH/used_parameters.yml ]; then
   echo "Cannot find used_parameters.yml in $BASE_PATH"
fi
parameters=($(parse_yaml $BASE_PATH/used_parameters.yml))

# Iterate over each string and find the parameters that we want
for i in "${parameters[@]}"
do
   if stringContain "Snapshots_basename=" "$i"; then
      SNAPSHOT_BASENAME=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_subdir=" "$i"; then
      SNAPSHOT_DIRECTORY=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_output_list=" "$i"; then
      OUTPUT_LIST=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_distributed=" "$i"; then
      SNAPSHOTS_ARE_DISTRIBUTED=`echo "$i" | cut -d'"' -f 2`
   fi
done

# Snapshots are, in principle, distributed. Each snapshot subfile is contained in its own directory,
# which has the same base name as the snapshot.
if [ $SNAPSHOTS_ARE_DISTRIBUTED == 1 ] ; then
   SNAPSHOT_BASEDIR=$SNAPSHOT_BASENAME
else
   $SNAPSHOT_BASEDIR=""
fi

# Check whether the snapshots are saved somewhere beyond the main directory (e.g. a 
# "snapshots" folder)
if [ $SNAPSHOT_DIRECTORY == "." ]; then
   SNAPSHOT_DIRECTORY=""

   # Check whether a snapshot folder exists, regardless of the default option. Added
   # because some users move the snapshots after being created, rather than using the
   # parameter file option
   if [ -d $BASE_PATH/snapshots ]; then
      SNAPSHOT_DIRECTORY="snapshots"
   fi
fi

# Snapshots that should be distributed are not, if one MPI rank was used. We defer this correction
# to be done here since we first needed to set SNAPSHOT_DIRECTORY
if [ $SNAPSHOTS_ARE_DISTRIBUTED == 1 ] ; then
   NUMBER_SUBDIR=$(find $BASE_PATH/$SNAPSHOT_DIRECTORY -maxdepth 1 -name "${SNAPSHOT_BASEDIR}_????" | wc -l)
   if [ $NUMBER_SUBDIR -eq 0 ]; then
      SNAPSHOT_BASEDIR=""
   fi
fi

# Get the number of snapshots we are expecting
NUM_OUTPUTS="$(grep -v '^#' $BASE_PATH/$OUTPUT_LIST | wc -l)"

# Check if the configuration file exists in the destination
if [ ! -f $HBT_FOLDER/config.txt ]; then

  echo "Copying configuration file to $HBT_FOLDER/config.txt"

  # Copy the template parameter file to the output folder
  cp ./template_config.txt $HBT_FOLDER/config.txt

  # Substitute the name of the simulation within the template
  sed -i "s@SIMULATION_PATH@${BASE_PATH}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the directory where snapshots are saved, if not in the
  # base simulation folder
  sed -i "s@SNAPSHOT_DIRECTORY@${SNAPSHOT_DIRECTORY}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the simulation within the template
  sed -i "s@OUTPUT_PATH@${HBT_FOLDER}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the snapshot
  sed -i "s@SNAPSHOT_BASEFILE@${SNAPSHOT_BASENAME}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the directory of snapshot subfiles, if applicable
  sed -i "s@SNAPSHOT_BASEDIR@${SNAPSHOT_BASEDIR}@g" $HBT_FOLDER/config.txt

  # Substitute the number of outputs within the HBT template
  sed -i "s@MAXIMUM_SNAPSHOT@$(($NUM_OUTPUTS - 1))@g" $HBT_FOLDER/config.txt

else
  echo "Configuration file already exists. Exiting to prevent accidental overwrites."
  exit 1
fi
