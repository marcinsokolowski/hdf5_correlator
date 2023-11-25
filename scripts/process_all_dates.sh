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

force=0
if [[ -n "$6" && "$6" != "-" ]]; then
   force=$6
fi

conversion_options=""
if [[ -n "$7" && "$7" != "-" ]]; then
   conversion_options="$7"   
fi

n_channels=32
if [[ -n "$8" && "$8" != "-" ]]; then
   n_channels=$8
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
                  processed_file=${dada_file%%dada}processed
            
                  if [[ -s $processed_file ]]; then
                     echo "File $dada_file already processed, in order to re-process remove file $processed_file"
                  else
                     echo "process_skalow_wide_bw_test.sh $dada_file $n_channels $channel 1 0 J0835-4510 $force \"$conversion_options\" > ${processed_file} 2>&1"
                     process_skalow_wide_bw_test.sh $dada_file $n_channels $channel 1 0 J0835-4510 $force "$conversion_options" > ${processed_file} 2>&1
                  fi
               done
               cd ..
            done
         fi
         
         echo "date > ${done_file}"
         date > ${done_file}
         
         cd ..
      done
   else
      echo "process_all_objects.sh \"${objects}\" ${conjugate}"
      process_all_objects.sh "${objects}" ${conjugate}   
   fi
      
   echo "date > ${done_file}"
   date > ${done_file}
#   fi
   cd ..
done

