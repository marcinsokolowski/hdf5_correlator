#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

freq_ch=204
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_ch=$2
fi


path=`which hdf5_to_dada_converter.py`

for datfile in `ls *.dat`
do
   outfile=${datfile%%.dat}_${object}.dada   
  
   echo "python $path ${datfile} --dat2dada --outfile=${outfile} --freq_ch=${freq_ch}"
   python $path ${datfile} --dat2dada --outfile=${outfile} --freq_ch=${freq_ch}
done
