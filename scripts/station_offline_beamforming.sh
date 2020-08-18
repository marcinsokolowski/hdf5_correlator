#!/bin/bash

station=eda2
if [[ -s /opt/aavs/config/station.yml ]]; then
   station=`awk -v station_section=0 '{if(index($1,":")>0 && NF==1){if(index($1,"station")>0 ){station_section=1;}else{station_section=0;}}if(station_section>0){if($1=="name:"){station_name=$2;gsub("\"","",station_name);gsub("\r","",station_name);print tolower(station_name);}}}' /opt/aavs/config/station.yml`
   echo "Station config file (or symbolik link) exists -> getting station_name = $station"
else
   echo "ERROR : /opt/aavs/config/station.yml file or symbolic link does not exist will use default station_name = $station or value passed in parameter -s"   
   exit;
fi


# WARNING : should be later in the script, but this value is required earlier :
pol_swap=0
if [[ $station == "eda2" ]]; then
   pol_swap=1
fi
if [[ -n "${10}" && "${10}" != "-" ]]; then
   pol_swap=${10}
fi

chan=204
if [[ -n "$1" && "$1" != "-" ]]; then
   chan=$1
fi

pol="XX"
pol_index=0
if [[ -n "$2" && "$2" != "-" ]]; then
   pol=$2
   if [[ $pol == "YY" ]]; then
      pol_index=1
   fi
fi

if [[ $pol_swap -gt 0 ]]; then
   echo "WARNING : swapping polarisations !!!"
   if [[ $pol_index -le 0 ]]; then
      # pol_index = 0 -> change to 1 
      pol_index=1
   else  
      # pol_index=1 -> change to 0 
      pol_index=0
   fi
fi

# trying to use last delays as default otherwise 256 zeros 
last_cal_path=/data/real_time_calibration/last_calibration/
last_cal_file=${last_cal_path}/chan_${chan}_selfcal_pha_${pol}.txt

delays=`echo 1 | awk '{for(i=0;i<256;i++){if(i<255){printf("0,");}else{printf("0");} }}'` # by default delays are zeros 
# delays can be overwritten by user-provided parameter
if [[ -n "$3" && "$3" != "-" ]]; then
   delays=$3
else
   if [[ -s ${last_cal_file} ]]; then
      delays=`~/aavs-calibration/station/calsol.sh ${last_cal_file}`
      echo "Last calibration file for channel $chan found : $last_cal_file -> delays are : $delays"
   fi
fi

if [[ -n "$4" && "$4" != "-" ]]; then
   station="$4"
fi

tag="lastcal_zenith_ch${chan}_${pol}"
if [[ -n "$5" && "$5" != "-" ]]; then
   tag="$5"
fi

merged_dir=`pwd`
if [[ -n "$6" && "$6" != "-" ]]; then
    merged_dir=$6
    cd $6/
fi


b=`basename $merged_dir`
if [[ $b != "merged" ]]; then
   echo "ERROR : off-line beamforming can only be executed on merged data (use option force to ignore directory name)"
   exit -1
fi

first_hdf5=`ls -tr *.hdf5 | head -1`
first_hdf5_date=`echo $first_hdf5 | cut -b 14-21`
dt=$first_hdf5_date
if [[ -n "$7" && "$7" != "-" ]]; then
   dt="$7"
fi

options=""
if [[ -n "$8" && "$8" != "-" ]]; then
   options=$8
fi

if [[ -n "$9" && "$9" != "-" ]]; then
   # passing specific path to calibration solutions to be used :
   last_cal_path=$9
   last_cal_file=${last_cal_path}/chan_${chan}_selfcal_pha_${pol}.txt
   
   echo "Using calibration solutions from file $last_cal_file"
   delays=`~/aavs-calibration/station/calsol.sh ${last_cal_file}`
   echo "Last calibration file for channel $chan found : $last_cal_file -> delays are : $delays"
fi

# pol_swap=1
# if [[ -n "${10}" && "${10}" != "-" ]]; then
#   pol_swap=${10}
# fi


echo "###################################################################"
echo "PARAMETERS:"
echo "###################################################################"
echo "station  = $station"
echo "pol_swap = $pol_swap"
echo "###################################################################"

# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

echo "nohup beamform_all.sh \"-C 1 -X $delays -q ${dt}_${station}_power_vs_time_ant -p $pol_index -u $pol ${options}\" ${dt}_${station}_256ant_${tag}.txt  > ${tag}.out 2>&1 &"
nohup beamform_all.sh "-C 1 -X $delays -q ${dt}_${station}_power_vs_time_ant -p $pol_index -u $pol ${options}" ${dt}_${station}_256ant_${tag}.txt  > ${tag}.out 2>&1 &
