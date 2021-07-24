#!/bin/bash

template="J* B*"
if [[ -n "$1" && "$1" != "-" ]]; then
   template="$1"
fi

for object in `ls -d ${template}`
do
   echo
   echo "----------------------------------------------------------------------------------------------------------------------------"
   date
   echo "dat2dada2_multifreq.sh $object - 0 > ${object}_processing.out 2>&1"
   dat2dada2_multifreq.sh $object - 0 > ${object}_processing.out 2>&1
done
