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

echo "######################################################################################"
echo "PARAMETERS of merge_all_dat.sh :"
echo "######################################################################################"
echo "do_process = $do_process"
echo "channel    = $channel"
echo "object     = $object"
echo "conjugate  = $conjugate"
echo "outdir     = $outdir"
echo "######################################################################################"


count=`ls *.dat | wc -l`
first_dat=`ls *.dat | head -1`
dat_path=`pwd`

if [[ $count -gt 1 ]]; then
   first_dat_b=${first_dat%%.dat}
   dadafile=${first_dat_b}.dada
   
   if [[ ! -s ${outdir}/${dadafile} ]]; then
      if [[ ! -s ${outdir}/${first_dat} ]]; then
         mkdir -p ${outdir}
         echo "cat *.dat > ${outdir}/${first_dat}"
         echo "It will take a bit of time ..."
         cat *.dat > ${outdir}/${first_dat} 
      else
         echo "INFO : file ${outdir}/${first_dat} already exists -> merging skipped"
      fi
   else
      echo "INFO : .dada file already exists -> no merging is required (remove .dada files to force re-mering)"
   fi
   
   if [[ $do_process -gt 0 ]]; then   
      if [[ -s channel.txt ]]; then
         channel_file=`cat channel.txt`
#         if [[ $channel_file != $channel ]]; then
#            echo "WARNING : channel in file channel.txt = $channel_file != $channel provided in the command line parameter -> please verify and re-start"
#            exit;
#         fi         
         channel=$channel_file
      else
         echo "WARNING : channel.txt does not exist using channel value = $channel from the command line, please verify that it is correct !"
         sleep 15
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
   if [[ $count -gt 0 ]]; then
      echo "DEBUG : just one .dat file found -> creating symbolic link"
      mkdir -p ${outdir}
      cd ${outdir}
      echo "ln -s ${dat_path}/${first_dat}"
      ln -s ${dat_path}/${first_dat}
      cd -
   else
      echo "WARNING : no .dat files in this directory -> nothing to be done ..."
      pwd
   fi
fi      
