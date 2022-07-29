#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

conjugate=1
if [[ -n "$3" && "$3" != "-" ]]; then
   conjugate=$3
fi

objects="J* B*"
if [[ -n "$4" && "$4" != "-" ]]; then
   objects="$4"
fi

multi_channel=0
if [[ -n "$5" && "$5" != "-" ]]; then
   multi_channel=$5
fi

force=1
if [[ -n "$6" && "$6" != "-" ]]; then
   force=$6
fi

conversion_options=""
if [[ -n "$7" && "$7" != "-" ]]; then
   conversion_options="$7"   
fi

n_channels=16
if [[ -n "$8" && "$8" != "-" ]]; then
   n_channels=$8
fi

dspsr_script=~/bin/dspsr_wrapper.sh
if [[ -n "$9" && "$9" != "-" ]]; then
   dspsr_script=$9
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt

echo "############################################"
echo "PARAMETERS:"
echo "############################################"
echo "datadir       = $datadir"
echo "template      = $template"
echo "conjugate     = $conjugate"
echo "objects       = $objects"
echo "multi_channel = $multi_channel"
echo "force         = $force"
echo "conversion_options = $conversion_options"
echo "n_channels    = $n_channels"
echo "dspsr_script  = $dspsr_script"
echo "############################################"
date

cd $datadir
pwd

for dir in `ls -d ${template}`
do
   echo
   echo "Processing $dir"
   pwd
   cd ${dir}
   pwd
   
#   if [[ -s ${done_file} ]]; then
#      echo "\t$dir already processed"
#   else
   if [[ $multi_channel -gt 0 ]]; then
      for subdir in `ls -d ${objects} 2>/dev/null`
      do
         object=`echo ${subdir} | cut -b 1-10`
         echo
         echo "INFO : processing subdirectory $subdir -> object = $object"
         
         if [[ -d ${subdir} ]]; then
            cd ${subdir}
         
            if [[ -s ${done_file} && $force -le 0 ]]; then
               echo "${subdir} / object ${object} already processed (remove file ${done_file} to repeat processing)"
            else         
               echo "INFO : processing ${subdir} / object ${object} ..."
               for channel in `ls -d ??? ?? 2>/dev/null`
               do
                  cd $channel
                  for dada_file in `ls *.dada 2>/dev/null`
                  do
                     # channel_1_1_1659060532.876270.dada
                     ch=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ch;}'`
                     channel_total=`echo "$channel $ch" | awk '{printf("%d\n",($1+$2));}'`
                     freq_mhz=`echo "$channel $ch" | awk '{printf("%.6f\n",($1+$2)*(400.00/512.00));}'`
                     ux=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ux;}'`
                     utc=`date -u -d "1970-01-01 UTC $ux seconds" +"%Y%m%dT%H%M%S"`
                     outfile=${utc}_ch${channel_total}
                     
                     echo ".dada file = $dada_file"
                     echo "ch = $channel + $ch = $channel_total -> freq = $freq_mhz [MHz]"
                     echo "ux = $ux -> utc  = $utc"
                     echo "outfile = $outfile"
                  
                     processed_file=${dada_file%%dada}processed
            
                     if [[ -s $processed_file && $force -le 0 ]]; then
                        echo "File $dada_file already processed, in order to re-process remove file $processed_file"
                     else
#                        echo "process_skalow_wide_bw_test.sh $dada_file $n_channels $channel 1 0 J0835-4510 $force \"$conversion_options\" > ${processed_file} 2>&1"
#                        process_skalow_wide_bw_test.sh $dada_file $n_channels $channel 1 0 J0835-4510 $force "$conversion_options" > ${processed_file} 2>&1
                         if [[ -s $dspsr_script ]]; then
                             echo "$dspsr_script $dada_file $object $channel"
                             $dspsr_script $dada_file $object $channel
                         else
                             if [[ ! -s ${object}.eph ]]; then
                                echo "psrcat -e ${object} > ${object}.eph"
                                psrcat -e ${object} > ${object}.eph
                             else
                                echo "Ephemeris file ${object}.eph already exists"
                             fi

                             size_mb=`du -smL ${dada_file} | awk '{print $1;}'`
                             echo "dspsr -E ${object}.eph -b 128 -F 256:D -U ${size_mb} -f ${freq_mhz} -O ${outfile} ${dada_file}"
                             dspsr -E ${object}.eph -b 128 -F 256:D -U ${size_mb} -f ${freq_mhz} -O ${outfile} ${dada_file}
                         fi
                     fi
                  done
                  cd ..
               done
            fi
         
            echo "date > ${done_file}"
            date > ${done_file}         
         cd ..
         else
            echo "WARNING : subdirectory $subdir does not exist -> skipped"
         fi
      done
   else
      # echo "process_all_objects.sh \"${objects}\" ${conjugate}"
      # process_all_objects.sh "${objects}" ${conjugate}   
      echo "ERROR : not expected option -> needs to be fixed"
   fi
      
   echo "date > ${done_file}"
   date > ${done_file}
#   fi
   cd ..
done

