#!/bin/bash

freq_channel=204
if [[ -n "$1" && "$1" != "-" ]]; then
   freq_channel=$1
fi
freq_mhz=`echo $freq_channel | awk '{printf("%.4f",$1*(400.00/512.00));}'`

station_name=eda2
if [[ -n "$2" && "$2" != "-" ]]; then
   station_name=$2
fi
station_name_lower=`echo $station_name | awk '{print tolower($1);}'`

tag=`date +%Y%m%d`
if [[ -n "$3" && "$3" != "-" ]]; then
   tag="$3"
fi

last_n_seconds=630720000
if [[ -n "$4" && "$4" != "-" ]]; then
   last_n_seconds=$4
fi

www_dir=aavs1-server:/exports/eda/${station_name_lower}/station_beam/
if [[ -n "$5" ]]; then
    www_dir=$5
fi

# AAVS_HOME = ~/Software on eda2-server
beam_scripts_path=$AAVS_HOME/hdf5_correlator/scripts/ # or on laptop ~/github/station_beam/processing/ or ~/Software on eda2 server 
if [[ -n "$6" && "$6" != "-" ]]; then
   beam_scripts_path=$6
fi

sleep_time=30
if [[ -n "$7" && "$7" != "-" ]]; then
   sleep_time=$7
fi

polarisation_swap=1 # in EDA2 (not in AAVS2) 
if [[ -n "$8" && "$8" != "-" ]]; then
   polarisation_swap=$8
fi


while [ 1 ];
do
   echo "~/Software/hdf5_correlator/scripts/process_station_beam.sh - ${freq_channel} ${station_name} ${tag} ${last_n_seconds} ${www_dir} ${beam_scripts_path} ${polarisation_swap}"
   ~/Software/hdf5_correlator/scripts/process_station_beam.sh - ${freq_channel} ${station_name} ${tag} ${last_n_seconds} ${www_dir} ${beam_scripts_path} ${polarisation_swap}
   
   sleep $sleep_time
   echo "sleep $sleep_time"
done

