#!/bin/bash

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

# Copies the catalogue sorting submission script into the HBT-HERONS folder, and creates
# the directories where the sorted catalogues and their logs will be saved.
function setup_sort_catalogues {
   local hbt_folder=$1
   cp ./submission_scripts/submit_sort_catalogues.sh $hbt_folder
   sed -i "s@CURRENT_PWD@${PWD}@g" $hbt_folder/submit_sort_catalogues.sh
   mkdir -p $hbt_folder/sorted_catalogues $hbt_folder/logs/sorted_catalogues
}
