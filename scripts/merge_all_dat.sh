#!/bin/bash

do_process=0
if [[ -n "$1" && "$1" != "-" ]]; then
   do_process=$1
fi

channel=410
if [[ -n "$2" && "$2" != "-" ]]; then
   channel=$2
fi

object=VELA
if [[ -n "$3" && "$3" != "-" ]]; then
   object=$3
fi

conjugate=1
if [[ -n "$4" && "$4" != "-" ]]; then
   conjugate=$4
fi

outdir=ALL/
if [[ -n "$5" && "$5" != "-" ]]; then
   outdir=$5
fi


count=`ls *.dat | wc -l`

if [[ $count -gt 0 ]]; then
   first_dat=`ls *.dat | head -1`

   if [[ ! -s ${outdir}/${first_dat} ]]; then
      mkdir -p ${outdir}
      echo "cat *.dat > ${outdir}/${first_dat}"
      echo "It will take a bit of time ..."
      cat *.dat > ${outdir}/${first_dat} 
   else
      echo "INFO : file ${outdir}/${first_dat} already exists -> merging skipped"
   fi
   
   if [[ $do_process -gt 0 ]]; then   
      if [[ -s channel.txt ]]; then
         channel=`cat channel.txt`
      fi      
   
      cd ${outdir}/
      first_dat_b=${first_dat%%.dat}
      echo "dat2dada2.sh ${object} ${channel} ${first_dat_b} 1 \"-F 256:D\" ${conjugate} - 1"
      sleep 5
      dat2dada2.sh ${object} ${channel} ${first_dat_b} 1 "-F 256:D" ${conjugate} - 1 # last 1 is to conjugate !   
      cd ..
   else
      echo "WARNING : full .dat file processing is not required, execute this line : merge_all_dat.sh 1"      
   fi
else
   echo "WARNING : no .dat files in this directory -> nothing to be done ..."
   pwd
fi      
