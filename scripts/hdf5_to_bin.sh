#!/bin/bash

hdf5_file_tile0=channel_cont_0_20190724_44084_0.hdf5
if [[ -n "$1" && "$1" != "-" ]]; then
   hdf5_file_tile0=$1
fi

convert_path=`which merge_n_hdf5_files.py`

echo "python $convert_path ${hdf5_file_tile0} 16"
python $convert_path ${hdf5_file_tile0} 16


merged_hdf5_file=`echo ${hdf5_file_tile0} | awk '{a=gsub("channel_cont_0_","channel_cont_",$1);print $1;}'`
echo "merged_hdf5_file = $merged_hdf5_file"

echo "hdf5dump.sh ${merged_hdf5_file}"
hdf5dump.sh ${merged_hdf5_file}
