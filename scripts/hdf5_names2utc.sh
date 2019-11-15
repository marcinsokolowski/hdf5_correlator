#!/bin/bash

hdf5_info_path=`which hdf5_info.py`
timestamp_from_hdf5file=0
inttime=1

for hdf5_file_tile0 in `ls *.hdf5`
do
    echo "---------------------------------------------------- $hdf5_file_tile0 ----------------------------------------------------"
    hdf5_info_file=${hdf5_file_tile0%%hdf5}hdf5_info

    echo "python $hdf5_info_path $hdf5_file_tile0 --inttime=${inttime} > ${hdf5_info_file}"
    python $hdf5_info_path ${hdf5_file_tile0} --inttime=${inttime} > ${hdf5_info_file}

    echo "---------------------------------------------"
    echo "hdf5 info :"
    echo "---------------------------------------------"
    cat ${hdf5_info_file}
    echo "---------------------------------------------"

    # convert seconds of the day in local time -> UTC and unixtime 
    if [[ $timestamp_from_hdf5file -gt 0 ]]; then
       dtm_ux=`grep "first timestamp" $hdf5_info_file | awk '{print $6;}'`
       dtm_local=`date -d "1970-01-01 UTC $dtm_ux seconds" +"%Y-%m-%d %T"`
       n_samples_per_hdf5_file=`grep n_integrations_per_uvfits $hdf5_info_file | awk '{print int($8);}'`

#       echo "INFO FROM hdf5_info.py output : dtm_ux = $dtm_ux (taken from $hdf5_info_file) -> dtm_local = $dtm_local , n_samples_per_hdf5_file = $n_samples_per_hdf5_file"
    else
       dtm_local=`echo $hdf5_file_tile0 | awk '{l=length($1);end=substr($1,l-20);dt=substr(end,1,8);daysec=substr(end,10,5);h=int(daysec/3600.00);m_dec=(((daysec/3600.00)-h)*60.00);m=int(m_dec);s=int((m_dec-int(m))*60.00);printf("%s %02d:%02d:%02d\n",dt,h,m,s);}'`
       # WORKS OK : date -u -d 'TZ="Australia/Perth" 20190724 10:25:07' +%s
       # dtm_string=`echo "'TZ=\"Australia/Perth\" ${dtm_local}'"`
       dtm_ux=`date -d "$dtm_local" +%s`
#       echo "dtm_local = $dtm_local -> $dtm_string"
    fi
    dtm_ut=`date -u -d "1970-01-01 UTC $dtm_ux seconds" +"%Y%m%d_%H%M%S"`
    lst=`ux2sid $dtm_ux mro | awk '{print $8;}'`
    
    echo "$hdf5_file_tile0 -> dtm_ut = $dtm_ut, ux = $dtm_ux , lst = $lst"
done
