#!/bin/bash

hdf5_file=`ls *.hdf5 | tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   hdf5_file=$1  
fi

antenna1=2 # was 0, but seems to be better to use 2 as it's also used in the normal calibration loop
if [[ -n "$2" && "$2" != "-" ]]; then
#    antenna1=$(($2*2))
   antenna1=$2
fi
antenna1_str=`echo $antenna1 | awk '{printf("%02d",$1);}'`

phase_norm=2
if [[ -n "$3" && "$3" != "-" ]]; then
   phase_norm=$3
fi
if [[ $phase_norm -le 0 ]]; then
    min_y=-360
    max_y=360
fi
if [[ $phase_norm -ge 2 ]]; then
    min_y=-180
    max_y=180
fi

n_ants=256
if [[ -n "$4" && "$4" != "-" ]]; then
   n_ants="$4"
fi

options=""
if [[ -n "$5" && "$5" != "-" ]]; then
   options="$5"
fi

antenna2=0

echo "######################################"
echo "PARAMETERS:"
echo "######################################"
echo "hdf5_file    = $hdf5_file"
echo "antenna1     = $antenna1"
echo "antenna2     = $antenna2"
echo "phase_norm   = $phase_norm ($min_y - $max_y)"
echo "n_ants       = $n_ants"
echo "options      = $options"
echo "######################################"
sleep 5

mkdir -p X/
cd X/
ln -s ../${hdf5_file}
echo "find_phases_from_corr! ${antenna1} 0 ${phase_norm} ${hdf5_file} ${n_ants} \"${options}\""
find_phases_from_corr! ${antenna1} 0 ${phase_norm} ${hdf5_file} ${n_ants} "${options}"
awk '{rest=substr($1,14);idx=index(rest,"_");ant2=substr(rest,5,3);print ant2" "$3;}' delta_phase.txt  > phase_vs_antenna_X.txt
# awk '{print substr($1,17,2)" "$3;}' delta_phase.txt  > phase_vs_antenna_X.txt
# echo corr_antenna_000_255.txt | awk '{rest=substr($1,14);idx=index(rest,"_");ant2=substr(rest,5,3);print rest" "idx" "ant2;}'
cd ../

mkdir -p Y/
cd Y/
ln -s ../${hdf5_file}
echo "find_phases_from_corr! ${antenna1} 1 ${phase_norm} ${hdf5_file} ${n_ants} \"${options}\""
find_phases_from_corr! ${antenna1} 1 ${phase_norm} ${hdf5_file} ${n_ants} "${options}"
awk '{rest=substr($1,14);idx=index(rest,"_");ant2=substr(rest,5,3);print ant2" "$3;}' delta_phase.txt  > phase_vs_antenna_Y.txt
# awk '{print substr($1,17,2)" "$3;}' delta_phase.txt  > phase_vs_antenna_Y.txt
cd ../

ln -fs X/phase_vs_antenna_X.txt
ln -fs Y/phase_vs_antenna_Y.txt

ls -tlr *.pkl

dt=`date +%Y%m%d`

echo "python $BIGHORNS/software/analysis/scripts/shell/eda/eda2/calibration.py --outfile=${dt}_calibration_coefficients.pkl"
python $BIGHORNS/software/analysis/scripts/shell/eda/eda2/calibration.py --outfile=${dt}_calibration_coefficients.pkl




