#!/bin/bash

dada_file=test.dada
if [[ -n "$1" && "$1" != "-" ]]; then
   dada_file=$1
fi

object=J0835-4510
if [[ -n "$2" && "$2" != "-" ]]; then
   object=$2
fi

channel=410
if [[ -n "$3" && "$3" != "-" ]]; then
   channel=$3
fi

ch=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ch;}'`
channel_total=`echo "$channel $ch" | awk '{printf("%d\n",($1+$2));}'`
freq_mhz=`echo "$channel $ch" | awk '{printf("%.6f\n",($1+$2)*(400.00/512.00));}'`
ux=`echo $dada_file | awk -F '_' '{ch=$2;ux=substr($4,1,17);print ux;}'`
utc=`date -u -d "1970-01-01 UTC $ux seconds" +"%Y%m%dT%H%M%S"`
outfile=${utc}_ch${channel_total}

echo ".dada file = $dada_file"
echo "ch = $channel + $ch = $channel_total -> freq = $freq_mhz [MHz]"
echo "ux = $ux -> utc  = $utc"
echo "outfile = $outfile"

if [[ ! -s ${object}.eph ]]; then
   echo "psrcat -e ${object} > ${object}.eph"
   psrcat -e ${object} > ${object}.eph
else
   echo "Ephemeris file ${object}.eph already exists"
fi

size_mb=`du -smL ${dada_file} | awk '{print $1;}'`

echo "dspsr -E ${object}.eph -b 128 -F 256:D -U ${size_mb} -f ${freq_mhz} -O ${outfile} ${dada_file}"
dspsr -E ${object}.eph -b 128 -F 256:D -U ${size_mb} -f ${freq_mhz} -O ${outfile} ${dada_file}


