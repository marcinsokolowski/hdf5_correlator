#!/bin/bash

# some checks :
psrcat_path=`which psrcat`
if [[ -n "$psrcat_path" ]]; then
   echo "OK : psrcat installed"
else
   echo "ERROR : psrcat not found -> cannot continue. Install from https://www.atnf.csiro.au/people/pulsar/psrcat/download.html"
   exit   
fi

channel_file="channel.txt"

for psrname in `ls -d J*`
do
   echo
   cd ${psrname}
   pwd
   date

   for ch in `ls -d ?? ???`
   do
      cd $ch   
      if [[ -s ${channel_file} ]]; then  
         channel_id=`cat ${channel_file}`
 
         echo "Processing pulsar $psrname observed at channel = $channel_id"

         # merge, but do not process :
         echo "merge_all_dat.sh 0"
         merge_all_dat.sh 0
         date
            
         cd ALL/
      
         # create eph file (if does not exist):
         if [[ -s ${psrname}.eph ]]; then
            echo "INFO : ephemeris file ${psrname}.eph already exists - nothing has to be done"       
         else
            echo "psrcat -e ${psrname} > ${psrname}.eph"
            psrcat -e ${psrname} > ${psrname}.eph
         fi
      
         echo "dat2dada2.sh ${psrname} ${channel_id} - 1 \"-F 256:D\" 1 - 1"
         dat2dada2.sh ${psrname} ${channel_id} - 1 "-F 256:D" 1 - 1 
      
      # ?      
#      pav -GCFDTp *.ar
 
      # The dynamic spectrum is
#      pav -GCDTpd *.ar

      # Or to fscrunch
#      pav -GTDTpd -f 16 *.ar
      
         if [[ -s ~/bin/pulsar_post_processing.sh ]]; then
            echo "File ~/bin/pulsar_post_processing.sh exists -> exacuting"
            chmod +x ~/bin/pulsar_post_processing.sh
         
            echo "~/bin/pulsar_post_processing.sh"
            ~/bin/pulsar_post_processing.sh
         else
            echo "WARNING : file ~/bin/pulsar_post_processing.sh does not exist -> no pulsar postprocessing performed on $psrname"
         fi
      
         cd ../
         date   
      else
         echo "ERROR : file $channel_file does not exist -> cannot process the pulsar"
      fi
      cd ..      
   done
   cd ..
done
