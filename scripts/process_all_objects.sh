#!/bin/bash

template="J* B*"
if [[ -n "$1" && "$1" != "-" ]]; then
   template="$1"
fi

conjugate=1
if [[ -n "$2" && "$2" != "-" ]]; then
   conjugate=$2
fi

for dir in `ls -d ${template}`
do
   object=`echo $dir | cut -b 1-10`    

   echo
   echo "----------------------------------------------------------------------------------------------------------------------------"
   date
   echo "dat2dada2_multifreq.sh $object - $conjugate >> ${object}_processing.out 2>&1"
   dat2dada2_multifreq.sh $object - $conjugate >> ${object}_processing.out 2>&1
done
