#!/bin/bash

do_process=0
if [[ -n "$1" && "$1" != "-" ]]; then
   do_process=$1
fi

count=`ls *.dat | wc -l`

if [[ $count -gt 0 ]]; then
   first_dat=`ls *.dat | head -1`

   mkdir -p ALL
   echo "cat *.dat ALL/${first_dat}"
   echo "It will take a bit of time ..."
   cat *.dat ALL/${first_dat}
   
   if [[ $do_process -gt 0 ]]; then   
      cd ALL/
      echo "dat2dada2.sh VELA 410 ${first_dat} 1 \"-F 256:D\" 1 - 1"
      sleep 5
      dat2dada2.sh VELA 410 ${first_dat} 1 "-F 256:D" 1 - 1 # last 1 is to conjugate !   
      cd ..
   else
      echo "WARNING : full .dat file processing is not required, execute this line : merge_all_dat.sh 1"      
   fi
else
   echo "WARNING : no .dat files in this directory -> nothing to be done ..."
   pwd
fi      
