#!/usr/bin/env bash
# shellcheck disable=SC1091

# set up variables
tfbs_beds=$1
vcfs_raw=$2
data_id_list=$3

# process each dataset ID
readarray -t data_ids < "$data_id_list"

for data in 



