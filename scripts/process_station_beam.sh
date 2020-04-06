#!/bin/bash

station_file=`ls *.hdf5  |tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   station_file=$1
fi
channel=4

tag=`date +%Y%m%d`
if [[ -n "$2" && "$2" != "-" ]]; then
   tag="$2"
fi

last_n_seconds=630720000
if [[ -n "$3" && "$3" != "-" ]]; then
   last_n_seconds=$3
fi

www_dir=aavs1-server:/exports/eda/eda2/station_beam/
if [[ -n "$4" ]]; then
    www_dir=$4
fi

freq_channel=204
if [[ -n "$5" && "$5" != "-" ]]; then
   freq_channel=$5
fi
freq_mhz=`echo $freq_channel | awk '{printf("%.4f",$1*(400.00/512.00));}'`

station_name=eda2
if [[ -n "$6" && "$6" != "-" ]]; then
   station_name=$6
fi


# AAVS_HOME = ~/Software on eda2-server
beam_scripts_path=$AAVS_HOME/hdf5_correlator/scripts/ # or on laptop ~/github/station_beam/processing/ or ~/Software on eda2 server 
if [[ -n "$7" && "$7" != "-" ]]; then
   beam_scripts_path=$7
fi



echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > x.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > x.out 2>&1

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > y.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} --freq_channel=${freq_channel} > y.out 2>&1

# temporary - TODO : add --freq_channel parameter to hdf5fits_station_beam.py
echo "cp ${tag}_power_vs_time_ch4_X.txt ${tag}_power_vs_time_ch${freq_channel}_X.txt"
cp ${tag}_power_vs_time_ch4_X.txt ${tag}_power_vs_time_ch${freq_channel}_X.txt
echo "cp ${tag}_power_vs_time_ch4_Y.txt ${tag}_power_vs_time_ch${freq_channel}_Y.txt"
cp ${tag}_power_vs_time_ch4_Y.txt ${tag}_power_vs_time_ch${freq_channel}_Y.txt

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

echo "python $beam_scripts_path/plot_power_vs_time.py ${tag}_power_vs_time_ch${freq_channel}"
python $beam_scripts_path/plot_power_vs_time.py ${tag}_power_vs_time_ch${freq_channel}


png_file=${tag}_power_vs_time_ch${freq_channel}.png

echo "cp images/${png_file} images/last_power_vs_time.png"
cp images/${png_file} images/last_power_vs_time.png

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


# root -l "plot_vs_time.C(\"power_vs_time_ch${channel}.txt\",-1e6,1e6,NULL,0,1000)"

