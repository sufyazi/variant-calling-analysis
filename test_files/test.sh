#!/usr/bin/env bash
# shellcheck disable=SC1091

which sem

arr=(1 2 3 4 5 6 7 8 9 10)
for i in "${arr[@]}"; do
    sem -j+0 echo "$i"
done