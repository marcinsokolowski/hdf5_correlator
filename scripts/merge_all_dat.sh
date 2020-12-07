#!/bin/bash

count=`ls *.dat | wc -l`

if [[ $count -gt 0 ]]; then
   first_dat=`ls *.dat | head -1`

   mkdir -p ALL
   cat *.dat ALL/${first_dat}
   cd ALL/
   echo "dat2dada2.sh VELA 410 ${first_dat} 1 \"-F 256:D\" 1 - 1"
   sleep 5
   dat2dada2.sh VELA 410 ${first_dat} 1 "-F 256:D" 1 - 1 # last 1 is to conjugate !   
   cd ..
else
   echo "WARNING : no .dat files in this directory -> nothing to be done ..."
   pwd
fi      
