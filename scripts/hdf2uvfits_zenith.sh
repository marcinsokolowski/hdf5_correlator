#!/bin/bash

# 1 channel : -b 1 
options=""
if [[ -n "$1" && "$1" != "-" ]]; then
   options="$1"
fi

hdf5_to_uvfits_all.sh -c -l -i 0.2831 -n 8140 -d "./" -N -z ${options}

