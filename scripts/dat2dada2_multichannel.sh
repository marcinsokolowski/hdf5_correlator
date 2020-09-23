#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

freq_ch=204
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_ch=$2
fi

do_dspsr=1
if [[ -n "$3" && "$3" != "-" ]]; then
   do_dspsr=$3
fi

dspsr_options=""
if [[ -n "$4" && "$4" != "-" ]]; then
   dspsr_options=$4
fi

conjugate=0
if [[ -n "$5" && "$5" != "-" ]]; then
   conjugate=$5
fi

force=0
if [[ -n "$6" && "$6" != "-" ]]; then
   force=$6
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
echo "###########################################################"



for ch in `ls -d ???`
do
   cd ${ch}/
   for dat_file in `ls *.dat`
   do  
      size_mb=`du -sm $dat_file`

      # only process files >500 MB 
      if [[ $size_mb -gt 500 ]];       
         echo "dat2dada2.sh ${object} ${freq_ch} ${dat_file} ${do_dspsr} "${dspsr_options}" $conjugate $force"
         dat2dada2.sh ${object} ${freq_ch} ${dat_file} ${do_dspsr} "${dspsr_options}" $conjugate $force
      else
         echo "WARNING : file $dat_file smaller than limit -> skipped"
      fi
   done
   cd ..
done

