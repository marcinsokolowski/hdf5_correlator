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

objects="J* B*"
if [[ -n "$4" && "$4" != "-" ]]; then
   objects="$4"
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt

echo "############################################"
echo "PARAMETERS:"
echo "############################################"
echo "datadir   = $datadir"
echo "template  = $template"
echo "conjugate = $conjugate"
echo "objects   = $objects"
echo "############################################"
date

cd $datadir
pwd

for dir in `ls -d ${template}`
do
   echo
   echo "Processing $dir"
   cd ${dir}
   pwd
   
#   if [[ -s ${done_file} ]]; then
#      echo "\t$dir already processed"
#   else
   echo "process_all_objects.sh \"${objects}\" ${conjugate}"
   process_all_objects.sh "${objects}" ${conjugate}   
      
   echo "date > ${done_file}"
   date > ${done_file}
#   fi
   cd -
done

