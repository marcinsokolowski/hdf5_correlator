#!/bin/bash

# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

# -H disables dumping of .bin files 
# -i 0.2831 is integration time in seconds when --channel_samples=262144  samples (1.08 usec each) are collected 
nohup bash `which hdf5_to_uvfits_all.sh` -i 0.2831 -n 8140 -d merged/ -H -S 0 > out 2>&1 &
