#!/bin/bash

hdf5_file=channel_cont_0_20190716_81307_0.hdf5
if [[ -n "$1" && "$1" != "-" ]]; then
   hdf5_file="$1"
fi

timesteps=-1
if [[ -n "$2" && "$2" != "-" ]]; then
   timesteps=$2
fi


# echo "hdf5_correlator ${hdf5_file} -a 0 -p 0 -n 2 -c 0  -L -x ${timesteps}"
# hdf5_correlator ${hdf5_file} -a 0 -p 0 -n 2 -c 0  -L -x ${timesteps}

binfile=${hdf5_file%%hdf5}bin

echo "hdf5_correlator ${hdf5_file} -d -o ${binfile}"
hdf5_correlator ${hdf5_file} -d -o ${binfile}

