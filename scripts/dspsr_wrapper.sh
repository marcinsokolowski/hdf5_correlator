#!/bin/bash

dada_file=test.dada
if [[ -n "$1" && "$1" != "-" ]]; then
   dada_file=$1
fi

object=J0835-4510
if [[ -n "$2" && "$2" != "-" ]]; then
   object=$2
fi

ch=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ch;}'`
freq_mhz=`echo $ch | awk '{printf("%.6f\n",$1*(400.00/512.00));}'`
ux=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ux;}'`
utc=`date -u -d "1970-01-01 UTC $1 seconds" +"%Y%m%dT%T"`
outfile=${utc}_ch${ch}.ar

echo ".dada file = $dada_file"
echo "ch = $ch -> freq = $freq_mhz [MHz]"
echo "ux = $ux -> utc  = $utc"
echo "outfile = $outfile"

if [[ ! -s ${object}.eph ]]; then
   echo "psrcat -e ${object} > ${object}.eph"
   psrcat -e ${object} > ${object}.eph
else
   echo "Ephemeris file ${object}.eph already exists"
fi

echo "dspsr -E ${object}.eph -b 64 -F 256:D -F 256:D -f ${freq_mhz} -O ${outfile} ${dada_file}"
dspsr -E ${object}.eph -b 64 -F 256:D -F 256:D -f ${freq_mhz} -O ${outfile} ${dada_file}


