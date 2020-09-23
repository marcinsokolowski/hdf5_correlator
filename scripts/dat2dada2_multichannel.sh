#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

do_dspsr=1
if [[ -n "$2" && "$2" != "-" ]]; then
   do_dspsr=$2
fi

dspsr_options=""
if [[ -n "$3" && "$3" != "-" ]]; then
   dspsr_options=$3
fi

conjugate=0
if [[ -n "$4" && "$4" != "-" ]]; then
   conjugate=$4
fi

force=0
if [[ -n "$5" && "$5" != "-" ]]; then
   force=$5
fi

auto=0
if [[ -n "$6" && "$6" != "-" ]]; then
   auto=$6
fi

echo "###########################################################"
echo "PARAMETERS:"
echo "###########################################################"
echo "object   = $object"
echo "freq_ch  = $freq_ch"
echo "prefix   = $prefix"
echo "do_dspsr = $do_dspsr"
echo "dspsr_options = $dspsr_options"
echo "conjugate = $conjugate"
echo "force    = $force"
echo "auto     = $auto"
echo "###########################################################"



for ch in `ls -d ???`
do
   cd ${ch}/
   ch_val=`echo $ch | awk '{printf("%d",$1);}'`
   
   for dat_file in `ls *.dat`
   do  
      size_mb=`du -sm $dat_file | awk '{print $1}'`
      b_dat_file=${dat_file%%.dat}

      # only process files >500 MB 
      if [[ $size_mb -gt 500 ]];       
      then
         echo "dat2dada2.sh ${object} ${ch_val} ${b_dat_file} ${do_dspsr} \"${dspsr_options}\" $conjugate $force $auto"
         dat2dada2.sh ${object} ${ch_val} ${b_dat_file} ${do_dspsr} "${dspsr_options}" $conjugate $force $auto
      else
         echo "WARNING : file $dat_file smaller than limit -> skipped"
      fi
   done
   cd ..
done

