#!/bin/bash

file=copied.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   file="$1"
fi

end_uxtime=$1

ux=`date +%s`
while [[ ! -s ${file} ]];
do
   echo "Waiting 10 seconds for file $file to be created ..."
   
   ux=`date +%s`
   sleep 10
done
