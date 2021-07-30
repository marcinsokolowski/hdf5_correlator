#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

conjugate=1
if [[ -n "$3" && "$3" != "-" ]]; then
   conjugate=$3
fi


echo "############################################"
echo "PARAMETERS:"
echo "############################################"
echo "datadir   = $datadir"
echo "template  = $template"
echo "conjugate = $conjugate"
echo "############################################"
date

cd $datadir
for dir in `ls -d ${template}`
do
   cd ${dir}
   echo "process_all_objects.sh \"J* B*\" ${conjugate}"
   process_all_objects.sh "J* B*" ${conjugate}   
   cd ..
done
