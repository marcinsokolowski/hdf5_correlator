#!/bin/bash

chan=204
if [[ -n "$1" && "$1" != "-" ]]; then
   chan=$1
fi

pol="XX"
if [[ -n "$2" && "$2" != "-" ]]; then
   pol=$2
fi

# trying to use last delays as default otherwise 256 zeros 
last_cal_path=/data/real_time_calibration/last_calibration/
last_cal_file=${last_cal_path}/chan_${chan}_selfcal_pha_${pol}.txt

delays=`echo 1 | awk '{for(i=0;i<256;i++){if(i<255){printf("0,");}else{printf("0");} }}'` # by default delays are zeros 
if [[ -s ${last_cal_file} ]]; then
   delays=`~/aavs-calibration/station/calsol.sh ${last_cal_file}`
   echo "Last calibration file for channel $chan found : $last_cal_file -> delays are : $delays"
fi

# delays can be overwritten by user-provided parameter
if [[ -n "$3" && "$3" != "-" ]]; then
   delays=$3
fi

station=eda2
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



# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

echo "nohup beamform_all.sh \"-C 1 -X $delays\" ${dt}_${station}_256ant_${tag}.txt  > ${tag}.out 2>&1 &"
nohup beamform_all.sh "-C 1 -X $delays" ${dt}_${station}_256ant_${tag}.txt  > ${tag}.out 2>&1 &
