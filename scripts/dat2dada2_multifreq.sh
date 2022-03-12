#!/bin/bash

name_template="J"
if [[ -n "$1" && "$1" != "-" ]]; then
   name_template="$1"
fi

phase_bins=64
if [[ -n "$2" && "$2" != "-" ]]; then
   phase_bins=$2
fi

conjugate=1
if [[ -n "$3" && "$3" != "-" ]]; then
   conjugate=$3
fi

outdir=ALL
if [[ $conjugate -le 0 ]]; then
   outdir=ALL_noConjugate
fi

done_file=processed.txt
force=0

echo "#######################################################"
echo "PARAMETERS:"
echo "#######################################################"
echo "name_template = $name_template"
echo "phase_bins    = $phase_bins"
echo "conjugate     = $conjugate"
echo "outdir        = $outdir"
echo "done_file     = $done_file"
echo "#######################################################"


# some checks :
echo "PSRCAT = $psrcat_path , PATH = $PATH"
psrcat_path=`which psrcat`
if [[ -n "$psrcat_path" ]]; then
   echo "OK : psrcat installed"
else
   echo "ERROR : psrcat not found -> cannot continue. Install from https://www.atnf.csiro.au/people/pulsar/psrcat/download.html"
   exit   
fi

channel_file="channel.txt"

for psrname in `ls -d ${name_template}*`
do
   echo
   cd ${psrname}
   pwd
   date
   
   echo "DEBUG : will process the following subdirectories (in this order) :"
   ls -d [0-9][0-9] [0-9][0-9]_* [0-9][0-9][0-9] [0-9][0-9][0-9]_* 2>/dev/null | sort -n -r

   # | sort -n -r : is to start from highest frequency channels (VELA best detected there):
   for ch in `ls -d [0-9][0-9] [0-9][0-9]_* [0-9][0-9][0-9] [0-9][0-9][0-9]_* 2>/dev/null | sort -n -r`
   do
      cd $ch   
      
      if [[ ! -s ${done_file} || $force -gt 0 ]]; then
         if [[ -s ${channel_file} ]]; then  
            channel_id=`cat ${channel_file}`
 
            echo "Processing pulsar $psrname observed at channel = $channel_id"

            # merge, but do not process :
            which merge_all_dat.sh
            echo "merge_all_dat.sh 0 - ${psrname} ${conjugate} ${outdir}"
            merge_all_dat.sh 0 - ${psrname} ${conjugate} ${outdir}
            date

            if [[ -d ${outdir} ]]; then            
               cd ${outdir}/
      
               # create eph file (if does not exist):
               if [[ -s ${psrname}.eph ]]; then
                  echo "INFO : ephemeris file ${psrname}.eph already exists - nothing has to be done"       
               else
                  echo "psrcat -e ${psrname} > ${psrname}.eph"
                  psrcat -e ${psrname} > ${psrname}.eph
               fi
      
               echo "dat2dada2.sh ${psrname} ${channel_id} - 1 \"-F 256:D\" ${conjugate} - 1 ${phase_bins}"
               dat2dada2.sh ${psrname} ${channel_id} - 1 "-F 256:D" ${conjugate} - 1 ${phase_bins}
      
               if [[ -s ~/bin/pulsar_post_processing.sh ]]; then
                  echo "File ~/bin/pulsar_post_processing.sh exists -> exacuting"
                  chmod +x ~/bin/pulsar_post_processing.sh
         
                  echo "~/bin/pulsar_post_processing.sh"
                  ~/bin/pulsar_post_processing.sh
               else
                  echo "WARNING : file ~/bin/pulsar_post_processing.sh does not exist -> no pulsar postprocessing performed on $psrname"
               fi
      
               cd ../
            
               echo "date > ${done_file}"
               date > ${done_file}
            else
               echo "ERROR : directory ${outdir} not created"
            fi
            date   
         else
            echo "ERROR : file $channel_file does not exist -> cannot process the pulsar"
         fi
      else
         echo "INFO : $done_file found -> data already processed (used force=1 to force re-processing"
      fi
      cd ..      
   done
   cd ..
done
