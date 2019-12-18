#!/bin/bash

ch=204
if [[ -n "$1" && "$1" != "-" ]]; then
   ch=$1   
fi

pol=XX
if [[ -n "$2" && "$2" != "-" ]]; then
   pol=$2
fi

dt="2???_??_??"
if [[ -n "$3" && "$3" != "-" ]]; then
  dt=$3
fi

caldir=/data/real_time_calibration/
if [[ -n "$4" && "$4" != "-" ]]; then
   caldir=$4
fi

# lastcal=`ls -d 2???_??_??-??:?? | tail -1` # I am assuming the SKA will not work beyond year 3000 ...

lastcal=`ls ${caldir}/${dt}-??:??/chan_${ch}_selfcal_pha_${pol}.txt | tail -1`
# lastcal_YY=`ls ${caldir}/2???_??_??-??:??/chan_${ch}_selfcal_pha_YY.txt | tail -1`

# use last calibration in the file :
tail -1 ${lastcal} | awk '{for(i=3;i<=NF;i++){printf("%.2f,",$i);}}' 




