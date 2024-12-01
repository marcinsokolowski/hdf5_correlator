#!/bin/bash

wait_for_file () {
  file=$1
  if [[ -n "$1" && "$1" != "-" ]]; then
     file="$1"
  fi

   ux=`date +%s`
   while [[ ! -s ${file} ]];
   do
      ux=`date +%s`
      echo "Waiting 10 seconds for file $file to be created (ux = $ux) ..."
      sleep 10
   done   
}


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

n_channels=40
if [[ -n "$4" && "$4" != "-" ]]; then
   n_channels=$4
fi

subdir_template="J?????????_flagants_ch${n_channels}* B???????_flagants_ch${n_channels}*"
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

for dir in `ls -trd ${template} | tail -5` # last 5 newest datasets only so that it does not get bogged down in old observations
do
   cd ${dir}
   pwd
   
   if [[ $force -gt 0 ]]; then      
      echo "force = $force -> removing file processed.txt"
      echo "rm -f processed.txt"
      rm -f processed.txt
   fi
   
   if [[ -s copied.txt ]]; then
      echo "INFO : data copied already -> can process"
   else
      echo "WARNING : data not fully copied yet -> waiting for it to be copied ..."
      echo "wait_for_file copied.txt"
      wait_for_file copied.txt      
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
