#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

freq_ch=204
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_ch=$2
fi

prefix="*"
if [[ -n "$3" && "$3" != "-" ]]; then
   prefix="$3"
fi

do_dspsr=1

path=`which hdf5_to_dada_converter.py`

for datfile in `ls ${prefix}.dat`
do
   unixtime=`echo $datfile | cut -b 11-25`
   echo "$datfile -> $unixtime - ok ?"
   sleep 2

   outfile=${datfile%%.dat}_${object}.dada   
   hdrfile=${datfile%%.dat}_${object}.hdr
  
#   echo "python $path ${datfile} --dat2dada --outfile=${outfile}"
#   python $path ${datfile} --dat2dada --outfile=${outfile}

   echo "python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch} --source=${object}"
   python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch} --source=${object}
   
   
   echo "cat ${hdrfile} ${datfile} > ${outfile}"
   cat ${hdrfile} ${datfile} > ${outfile}
   
   if [[ $do_dspsr -gt 0 ]]; then
      echo "dspsr -E 1752.eph -b 64 -U 600 ${outfile}"
      dspsr -E 1752.eph -b 64 -U 600 ${outfile}
   
      last_ar=`ls -tr *.ar | tail -1`
   
      echo "psrplot -p flux -D /xs $last_ar"
      psrplot -p flux -D /xs $last_ar
   else
      echo "WARNING : dspsr is not required"
   fi
done
