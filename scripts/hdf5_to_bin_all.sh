#!/bin/bash

do_correlation=0
bin2lfiles=0
ra_hrs=
dec_degs=
lfile_converter_path="./"

function print_usage {
  echo "Script merges all hdf5 files from 16 tiles to a single one, converts them into bin files and depending on options runs correlation and lfiles to uvfits conversion"
  echo "Usage: "
  echo "hdf5_to_bin_all.sh [options]"
  echo "    -c Enables correlation Default: $do_correlation"
  echo "    -l Converts l-files to uvfits files Default: $bin2lfiles"
  echo "    -R ra_hours   Default: use zenith"
  echo "    -D dec_degs   Default: use zenith"
  exit
}


# parse command-line args
if [ $# -lt 1 ] ; then print_usage ; fi
while getopts "hclR:D:" opt; do
  case $opt in
    h)
        print_usage
        ;;
    c)
        do_correlation=1
        ;;
    l)
        bin2lfiles=1
        ;;
    R)
        ra_hrs=$OPTARG
        useradec=1
        ;;
    D)
        dec_degs=$OPTARG
        ;;
    \?)
      echo "Invalid option: -$OPTARG" 1>&2
      print_usage
      ;;
  esac
done
shift $(expr $OPTIND - 1 )


for hdf5_tile0_file in `ls channel_cont_0_????????_*_0.hdf5`
do
    merged_hdf5_file=`echo ${hdf5_file_tile0} | awk '{a=gsub("channel_cont_0_","channel_cont_",$1);print $1;}'`
    merged_bin_file=${merged_hdf5_file%%hdf5}bin

    echo "hdf5_to_bin.sh $hdf5_tile0_file"
    hdf5_to_bin.sh $hdf5_tile0_file       
    
    dt=
    
    lfile_base=


    if [[ $do_correlation -gt 0 ]]; then
#     
#       /home/rwayth/bin/corr_gpu_complex2 -c 32 -n 512 -a 2894 -i 20190724_041444_eda2_ch32_ant256_midday_avg2894.bin -o 20190724_041444_eda2_ch32_ant256_midday_avg2894 -w 10
    fi 

    if [[ $bin2lfiles -gt 0 ]]; then
       radec_string=""
       if [[ -n "$ra_hrs" && -n "$dec_degs" ]]; then
          radec_string="-R $ra_hrs -D $dec_degs"
          echo "radec_string = $radec_string"
       fi
       echo "${lfile_converter_path}/Lfile2uvfits_eda.sh $lfile_base $radec_string"
       ${lfile_converter_path}/Lfile2uvfits_eda.sh $lfile_base $radec_string
    else 
       echo "WARNING : conversion from .bin -> lfiles is not required"
    fi    
done

