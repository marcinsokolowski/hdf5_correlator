#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

freq_start=410
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_start=$3
fi

objects="J* B*"
if [[ -n "$4" && "$4" != "-" ]]; then
   objects="$4"
fi

multi_channel=1
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

n_channels=128
if [[ -n "$8" && "$8" != "-" ]]; then
   n_channels=$8
fi

period=0.0894189988
if [[ -n "$9" && "$9" != "-" ]]; then
   period=$9
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt

echo "############################################"
echo "PARAMETERS:"
echo "############################################"
echo "datadir       = $datadir"
echo "template      = $template"
echo "freq_start    = $freq_start"
echo "objects       = $objects"
echo "multi_channel = $multi_channel"
echo "force         = $force"
echo "conversion_options = $conversion_options"
echo "n_channels    = $n_channels"
echo "period        = $period"
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
            cd ${subdir}/${freq_start}/
         
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
                         # skalow_spectrometer test.dada -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F 410 -N 128 -O dynspec -a 7 -P 0.0894189988 -D 2
                         echo "skalow_spectrometer $dada_file -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F ${channel} -N ${n_channels} -O dynspec -a 7 -P ${period} -D 2"
                         skalow_spectrometer $dada_file -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F ${channel} -N ${n_channels} -O dynspec -a 7 -P ${period} -D 2
                     fi
                  done
                  cd ..
               done
            fi
         
            echo "date > ${done_file}"
            date > ${done_file}         
            cd ../..
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

