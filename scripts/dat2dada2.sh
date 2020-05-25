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
if [[ -n "$4" && "$4" != "-" ]]; then
   do_dspsr=$4
fi

dspsr_options=""
if [[ -n "$5" && "$5" != "-" ]]; then
   dspsr_options=$5
fi

force=0

eph_dir=~/github/hdf5_correlator/scripts/config/dspsr/

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

   if [[ ! -s ${outfile} || $force -gt 0 ]]; then
      echo "python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch} --source=${object}"
      python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch} --source=${object}
   
      size_mb=`du -sm ${datfile} | awk '{print $1;}'`
      echo "size_mb = $size_mb"
   
      echo "cat ${hdrfile} ${datfile} > ${outfile}"
      cat ${hdrfile} ${datfile} > ${outfile}
   
      if [[ $do_dspsr -gt 0 ]]; then
         if [[ ! -s ${object}.eph ]]; then
            echo "cp ${eph_dir}/${object}.eph ."
            cp ${eph_dir}/${object}.eph .
         fi
   
         if [[ -s ${object}.eph ]]; then
            echo "dspsr -E ${object}.eph -b 64 -U $size_mb ${dspsr_options} ${outfile}"
            dspsr -E ${object}.eph -b 64 -U $size_mb ${dspsr_options} ${outfile}
   
            last_ar=`ls -tr *.ar | tail -1`
   
            echo "psrplot -p flux -D /xs $last_ar"
            psrplot -p flux -D /xs $last_ar
            
            echo "pav -G -DTp -N1,1 2 $last_ar"
            pav -G -DTp -N1,1 2 $last_ar
         else
            echo "WARNING : missing file ${object}.eph , cannot find local version neither in ${eph_dir} - please fix it and re-run dspsr"
            echo "dspsr -E ${object}.eph -b 64 -U 600 ${dspsr_options} ${outfile}"
            echo "and : psrplot -p flux -D /xs $last_ar"
         fi
      else
         echo "WARNING : dspsr is not required"
      fi
   else
      echo "WARNING : dada file ${outfile} already exists -> skipped (enable force=1 in order to re-process)"
   fi
done

