#!/bin/bash

channel=0
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

n_channels=2
if [[ -n "$2" && "$2" != "-" ]]; then
   n_channels=$2
fi


mkdir -p merged_channel${channel}/

path=`which merge_hdf5_n_chan.py`
if [[ -n "$path" ]]; then
   echo "OK : script merge_hdf5_n_chan.py found in $path"
else
   echo "ERROR : script merge_hdf5_n_chan.py not found on PATH, trying doing this:"
   echo "export PATH=~/github/hdf5_correlator/scripts:\$PATH"
   echo "or"
   echo "export PATH=~/Software/hdf5_correlator/scripts:\$PATH"
   echo "or both, or wherever hdf5_correlator package is located"
fi

for hdf5file in `ls channel_cont_0_*.hdf5`
do
   echo "python $path $hdf5file 16 - --outdir=merged_channel${channel}/ --out_channel=${channel} --n_channels=${n_channels}"
   python $path $hdf5file 16 - --outdir=merged_channel${channel}/ --out_channel=${channel} --n_channels=${n_channels}
done
