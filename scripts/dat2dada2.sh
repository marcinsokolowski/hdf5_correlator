#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

path=`which hdf5_to_dada_converter.py`

for datfile in `ls *.dat`
do
   outfile=${datfile%%.dat}_${object}.dada   
   hdrfile=${datfile%%.dat}_${object}.hdr
  
#   echo "python $path ${datfile} --dat2dada --outfile=${outfile}"
#   python $path ${datfile} --dat2dada --outfile=${outfile}

   echo "python $path ${datfile} --psrdadahdrs --outfile=${hdrfile}"
   python $path ${datfile} --psrdadahdrs --outfile=${hdrfile}
   
   
   echo "cat ${hdrfile} ${datfile} > ${outfile}"
   cat ${hdrfile} ${datfile} > ${outfile}
done
