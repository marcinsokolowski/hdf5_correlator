#!/bin/bash

channel=-1
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

if [[ $channel -lt 0 ]]; then
   echo "ERROR : channel not specified ! Try again as : convert.sh CHANNEL"
   exit
fi

# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

# -H disables dumping of .bin files 
# -i 0.2831 is integration time in seconds when --channel_samples=262144  samples (1.08 usec each) are collected 
pwd
bash `which hdf5_to_uvfits_all.sh` -i 0.2831 -n 8140 -d merged/ -H -S 0
date

cd merged/
pwd
echo "/usr/local/bin/hdf2uvfits_sun_aavs2.sh \"-f ${channel}\""
/usr/local/bin/hdf2uvfits_sun_aavs2.sh "-f ${channel}"

date




