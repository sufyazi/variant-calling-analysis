#!/bin/bash

# assign the first argument to a variable
json_object_file=$1

# Initialize an associative array
declare -A id_map

# Loop through each JSON object from the pipe
while read -r json; do
  submitter_id=$(echo "$json" | jq -r '.submitter_id')
  aliquot_id=$(echo "$json" | jq -r '.aliquot_id')

  if [[ $submitter_id == *aliquot ]]  ; then
    echo "$submitter_id contains aliquot"
    # Populate the associative array
    id_map["$submitter_id"]=$aliquot_id
  fi

done < <(cat "$json_object_file")

# Print the associative array
for key in "${!id_map[@]}"; do
  echo "Submitter ID: $key Aliquot ID: ${id_map[$key]}"
done