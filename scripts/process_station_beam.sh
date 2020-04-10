#!/bin/bash

station_file=`ls *.hdf5  |tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   station_file=$1
fi
channel=4

# temporary solution - as it's not possible to share HDF5 file ???
# create temporary copy of station HDF5 file :
echo "cp ${station_file} station_beam_temporary.hdf5"
cp ${station_file} station_beam_temporary.hdf5
station_file=station_beam_temporary.hdf5

freq_channel=204
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_channel=$2
fi
freq_mhz=`echo $freq_channel | awk '{printf("%.4f",$1*(400.00/512.00));}'`

station_name=eda2
if [[ -n "$3" && "$3" != "-" ]]; then
   station_name=$3
fi

tag=`date +%Y%m%d`
if [[ -n "$4" && "$4" != "-" ]]; then
   tag="$4"
fi

last_n_seconds=630720000
if [[ -n "$5" && "$5" != "-" ]]; then
   last_n_seconds=$5
fi

www_dir=aavs1-server:/exports/eda/eda2/station_beam/
if [[ -n "$6" ]]; then
    www_dir=$6
fi

# AAVS_HOME = ~/Software on eda2-server
beam_scripts_path=$AAVS_HOME/hdf5_correlator/scripts/ # or on laptop ~/github/station_beam/processing/ or ~/Software on eda2 server 
if [[ -n "$7" && "$7" != "-" ]]; then
   beam_scripts_path=$7
fi

polarisation_swap=1 # in EDA2 (not in AAVS2) 
if [[ -n "$8" && "$8" != "-" ]]; then
   polarisation_swap=$8
fi

echo "ls -al station*.hdf5"
ls -al station*.hdf5

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > x.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > x.out 2>&1

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > y.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > y.out 2>&1

echo "ls -al ${tag}_power_vs_time_ch*"
ls -al ${tag}_power_vs_time_ch*

# Polarisation swap :
if [[ $polarisation_swap -gt 0 ]]; then
   echo "mv ${tag}_power_vs_time_ch${freq_channel}_X.txt ${tag}_power_vs_time_ch${freq_channel}_Y.tmp"
   mv ${tag}_power_vs_time_ch${freq_channel}_X.txt ${tag}_power_vs_time_ch${freq_channel}_Y.tmp

   echo "mv ${tag}_power_vs_time_ch${freq_channel}_Y.txt ${tag}_power_vs_time_ch${freq_channel}_X.txt"
   mv ${tag}_power_vs_time_ch${freq_channel}_Y.txt ${tag}_power_vs_time_ch${freq_channel}_X.txt

   echo "mv ${tag}_power_vs_time_ch${freq_channel}_Y.tmp ${tag}_power_vs_time_ch${freq_channel}_Y.txt"
   mv ${tag}_power_vs_time_ch${freq_channel}_Y.tmp ${tag}_power_vs_time_ch${freq_channel}_Y.txt
else
   echo "INFO : no polarisation swap required"
fi

ls ${tag}_power_vs_time_ch*.txt > list
first_file=`head --lines=1 list`

mkdir -p images/

root_path=`which root`
if [[ -n $root_path ]]; then
   echo "Plotting in ROOT ..."
   root -l "plotNfiles_vs_time_scaling.C(\"list\",-1,1,-1,1.00)"
else
   echo "WARNING : ROOT cern software not installed, skipping plot in ROOT (no worries !)"
fi

# When no --y_min and --y_max specified -> AUTO-SCALE 
# was --y_min=0 --y_max=5000 or --y_min=0 --y_max=1000
echo "python $beam_scripts_path/plot_power_vs_time.py ${tag}_power_vs_time_ch${freq_channel} --comment=\"${tag}_power_vs_time_ch${freq_channel} (pol. swapped)\""
python $beam_scripts_path/plot_power_vs_time.py ${tag}_power_vs_time_ch${freq_channel} --comment="${tag}_power_vs_time_ch${freq_channel} (pol. swapped)"


png_file=${tag}_power_vs_time_ch${freq_channel}.png

echo "cp images/${png_file} images/last_power_vs_time.png"
cp images/${png_file} images/last_power_vs_time.png

dtm=`date +%Y%m%d%M%S`
echo "cp images/${png_file} images/${dtm}.png"
cp images/${png_file} images/${dtm}.png

if [[ -n ${www_dir} && ${www_dir} != "-" ]]; then
   echo "INFO : Copying image to ${www_dir}/"

   echo "rsync -avP images/${png_file} ${www_dir}/"
   rsync -avP images/${png_file} ${www_dir}/

   echo "rsync -avP images/last_power_vs_time.png ${www_dir}/"
   rsync -avP images/last_power_vs_time.png ${www_dir}/   
   
   echo "cp $beam_scripts_path/html/power.template power.html"
   cp $beam_scripts_path/html/power.template power.html
   
   station_name_uppper=`echo $station_name | awk '{print toupper($1);}'`     
   
   echo "sed -i 's/FREQ_CHANNEL/${freq_channel}/g' power.html" > sed!
   echo "sed -i 's/FREQ_VALUE_MHZ/${freq_mhz}/g' power.html" >> sed!
   echo "sed -i 's/STATION_NAME/${station_name_uppper}/g' power.html" >> sed!
   chmod +x sed!
   ./sed!

   echo "rsync -avP power.html ${www_dir}/"
   rsync -avP power.html ${www_dir}/
else
   echo "WARNING : www_dir not specified, image not copied anywhere"
fi   

# remove temporary copy 
echo "rm -f station_beam_temporary.hdf5"
rm -f station_beam_temporary.hdf5

# root -l "plot_vs_time.C(\"power_vs_time_ch${channel}.txt\",-1e6,1e6,NULL,0,1000)"

