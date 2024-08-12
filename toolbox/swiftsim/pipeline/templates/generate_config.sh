#!/bin/bash
set -e

# Parses a YAML file and returns each option as its own string. The levels of each
# option are reflected by concatenating underscores between names.
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Identifies whether a string contains a substring.
stringContain() { case $2 in *$1* ) return 0;; *) return 1;; esac ;}

# This is where the outputs should be placed
BASE_PATH=$1
HBT_FOLDER=$2

# Parse the SWIFT config file
if [ ! -f $BASE_PATH/used_parameters.yml ]; then
   echo "Cannot find used_parameters.yml in $1"
fi
parameters=($(parse_yaml $BASE_PATH/used_parameters.yml))

# Iterate over each string and find the parameters that we want
for i in "${parameters[@]}"
do
   if stringContain "Snapshots_basename=" "$i"; then
      SNAPSHOT_BASENAME=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_subdir=" "$i"; then
      SNAPSHOT_BASEDIR=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_output_list=" "$i"; then
      OUTPUT_LIST=`echo "$i" | cut -d'"' -f 2`
   fi

   if stringContain "Snapshots_distributed=" "$i"; then
      SNAPSHOTS_ARE_DISTRIBUTED=`echo "$i" | cut -d'"' -f 2`
   fi
done

if [ $SNAPSHOTS_ARE_DISTRIBUTED == 1 ] ; then
   # Snapshots are distributed, but a directory base name has not been specified. It
   # defaults to the same base name as the snapshots.
   if [ $SNAPSHOT_BASEDIR == "." ]; then
      echo "SNAPSHOTS ARE DISTRIBUTED but we are using default"
      SNAPSHOT_BASEDIR=$SNAPSHOT_BASENAME
   fi
   # Snapshots are not distributed, so we should leave BaseDir blank
else
   $SNAPSHOT_BASEDIR=""
fi
echo $SNAPSHOT_BASENAME $OUTPUT_LIST $SNAPSHOTS_ARE_DISTRIBUTED $SNAPSHOT_BASEDIR

# Get the number of snapshots we are expecting
NUM_OUTPUTS="$(grep -v '^#' $BASE_PATH/$OUTPUT_LIST | wc -l)"

# Check if the configuration file exists in the destination
if [ ! -f $HBT_FOLDER/config.txt ]; then

  echo "Copying configuration file to $HBT_FOLDER/config.txt"

  # Copy the template parameter file to the output folder
  cp ./template_config.txt $HBT_FOLDER/config.txt

  # Substitute the name of the simulation within the template
  sed -i "s@SIMULATION_PATH@${BASE_PATH}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the simulation within the template
  sed -i "s@OUTPUT_PATH@${HBT_FOLDER}@g" $HBT_FOLDER/config.txt

  # Substitute the name of the subdirectory, if any
  sed -i "s@SNAPSHOT_BASEFILE@${SNAPSHOT_BASENAME}@g" $HBT_FOLDER/config.txt
  #
  # Substitute the name of the subdirectory, if any
  sed -i "s@SNAPSHOT_BASEDIR@${SNAPSHOT_BASEDIR}@g" $HBT_FOLDER/config.txt

  # Substitute the number of outputs within the HBT template
  sed -i "s@MAXIMUM_SNAPSHOT@$(($NUM_OUTPUTS - 1))@g" $HBT_FOLDER/config.txt

else
  echo "Configuration file already exists. Exiting to prevent accidental overwrites."
  exit 1
fi