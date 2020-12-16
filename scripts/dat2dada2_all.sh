#!/bin/bash

channel_file="channel.txt"

for psrname in `ls -d J* V*`
do
   echo
   cd ${psrname}
   pwd
   date
   
   if [[ -s ${channel_file} ]]; then  
      channel_id=`cat ${channel_file}`
 
      echo "Processing pulsar $psrname observed at channel = $channel_id"

      # merge, but do not process :
      echo "merge_all_dat.sh 0"
      merge_all_dat.sh 0
      date
            
      cd ALL/
      echo "dat2dada2.sh ${psrname} ${channel_id} - 1 \"-F 256:D\" 1 - 1"
      dat2dada2.sh ${psrname} ${channel_id} - 1 "-F 256:D" 1 - 1 
      
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
