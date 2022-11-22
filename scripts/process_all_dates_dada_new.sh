#!/bin/bash

datadir="/data_archive/"
if [[ -n "$1" && "$1" != "-" ]]; then
   datdir="$1"
fi

template="202[2-3]_??_??_pulsars_vela"
if [[ -n "$2" && "$2" != "-" ]]; then
   template="$2"
fi

force=0
if [[ -n "$3" && "$3" != "-" ]]; then
   force=$3
fi

echo "##############################################################"
echo "PARAMETERS :"
echo "##############################################################"
echo "datadir  = $datadir"
echo "template = $template"
echo "force    = $force"
echo "##############################################################"


export PATH=/home/msok/github/hdf5_correlator/scripts/:$PATH

cd ${datadir}

for dir in `ls -d ${template}`
do
   cd ${dir}
   
   if [[ $force -gt 0 ]]; then      
      echo "force = $force -> removing file processed.txt"
      echo "rm -f processed.txt"
      rm -f processed.txt
   fi
   
   if [[ -s processed.txt ]]; then         
      echo "INFO : data in ${dir} already processed -> skipped"
   else
      echo "process_all_dates_dada.sh /data_archive/ ${dir} 1 "J?????????_flagants_16ch*" 1 1 - 16 > process.out 2>&1"
      process_all_dates_dada.sh /data_archive/ ${dir} 1 "J?????????_flagants_16ch*" 1 1 - 16 > process.out 2>&1 
   
      date >> processed.txt
      echo "Processing of $dir done at :"
      date 
   fi
   cd ..
done
