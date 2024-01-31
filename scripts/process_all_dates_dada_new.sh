#!/bin/bash

datadir="/data_archive/"
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir="$1"
fi

template="202[2-3]_??_??_pulsars_vela"
if [[ -n "$2" && "$2" != "-" ]]; then
   template="$2"
fi

force=0
if [[ -n "$3" && "$3" != "-" ]]; then
   force=$3
fi

n_channels=16
if [[ -n "$4" && "$4" != "-" ]]; then
   n_channels=$4
fi

subdir_template="J?????????_flagants_${n_channels}ch* B???????_flagants_${n_channels}ch*"
if [[ -n "$5" && "$5" != "-" ]]; then
   subdir_template="$5"
fi


echo "##############################################################"
echo "PARAMETERS :"
echo "##############################################################"
echo "datadir    = $datadir"
echo "template   = $template"
echo "force      = $force"
echo "n_channels = $n_channels"
echo "subdir_template = $subdir_template"
echo "##############################################################"


# PULSARS :
export PATH=/home/msok/github/hdf5_correlator/scripts/:/usr/local/bin/:/home/aavs/bin/psrcat_tar:$PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/data1/msok/anaconda2/lib/:$LD_LIBRARY_PATH # for gfortran in anaconda2
export TEMPO=/home/aavs/Software/pulsars/tempo
export TEMPO2=/usr/share/tempo2/


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
      echo "Processing started at :" > processed.txt
      date >> processed.txt
      echo "process_all_dates_dada.sh ${datadir} ${dir} 1 \"${subdir_template}\" 1 1 - ${n_channels} > process.out 2>&1"
      process_all_dates_dada.sh ${datadir} ${dir} 1 "${subdir_template}" 1 1 - ${n_channels} > process.out 2>&1 
   
      echo "Processing finished at :" > processed.txt
      date >> processed.txt
      echo "Processing of $dir done at :"
      date 
   fi
   cd ..
done
