#!/bin/bash

station_name=eda2

if [[ -s /opt/aavs/config/station.yml ]]; then
   station_name=`awk -v station_section=0 '{if(index($1,":")>0 && NF==1){if(index($1,"station")>0 ){station_section=1;}else{station_section=0;}}if(station_section>0){if($1=="name:"){station_name=$2;gsub("\"","",station_name);print tolower(station_name);}}}' /opt/aavs/config/station.yml`   
#   echo "Station config file (or symbolik link) exists -> getting station_name = $station_name"
fi

echo $station_name
