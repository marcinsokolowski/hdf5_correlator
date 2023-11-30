#!/bin/bash

dat_file=`ls *.dat | tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   dat_file=$1
fi
# ext=`echo $dat_file | cut -d'.' -f2`
ext=`echo $dat_file | awk '{ext=$1;i=index(ext,".");while(i>0){ext_new=substr(ext,i+1);ext=ext_new;i=index(ext,".");}}END{print ext;}'`

n_channels=32
if [[ -n "$2" && "$2" != "-" ]]; then
   n_channels=$2
fi

ch_start=112
if [[ -n "$3" && "$3" != "-" ]]; then
   ch_start=$3
fi

process_channel=0
if [[ -n "$4" && "$4" != "-" ]]; then
   process_channel=$4
fi

test_pattern_generator=0
if [[ -n "$5" && "$5" != "-" ]]; then
   test_pattern_generator=$5
fi

object="J0835-4510"
if [[ -n "$6" && "$6" != "-" ]]; then
   object="$6"
fi

force=0
if [[ -n "$7" && "$7" != "-" ]]; then
   force=$7
fi

conversion_options=""
if [[ -n "$8" && "$8" != "-" ]]; then
   conversion_options="$8"   
fi

do_psr_processing=1
if [[ -n "$9" && "$9" != "-" ]]; then
   do_psr_processing=$9
fi


start_byte=0
if [[ $ext == "dada" ]]; then
   # to skip .dada header 
   start_byte=4096

   echo "head --bytes=${start_byte} $dat_file > dada_header.txt"   
   head --bytes=${start_byte} $dat_file > dada_header.txt
   
   if [[ ! -s header_template.txt ]]; then
      echo "cp dada_header.txt header_template.txt"
      cp dada_header.txt header_template.txt
   else
      echo "INFO : header_template.txt file already exists"
   fi
fi

echo "###############################################################"
echo "PARAMETERS :"
echo "###############################################################"
echo "dat_file   = $dat_file (extension = $ext -> start_byte = $start_byte)"
echo "n_channels = $n_channels"
echo "ch_start   = $ch_start"
echo "process_channel = $process_channel"
echo "test_pattern_generator = $test_pattern_generator"
echo "force      = $force"
echo "proces     = $process_channel"
echo "conversion_options = $conversion_options"
echo "do_psr_processing = $do_psr_processing"
echo "###############################################################"

path_dada_converter=`which dada_header_converter.py`

ch=0
while [[ $ch -lt ${n_channels} ]];
do
   ch_full=$(($ch_start+$ch))
   ch_str=`echo $ch_full | awk '{printf("%03d",$1);}'`
   
   cnt_ar=`ls ${ch_str}/*.ar 2>/dev/null | wc -l`
   echo "Number of ar files = $cnt_ar"
   
   if [[ ! -s ${ch_str}/test.dada || $cnt_ar -le 0 || $force -gt 0 ]]; then
      mkdir -p ${ch_str}

      if [[ -s ${ch_str}/test.dat ]]; then
         echo "INFO : file ${ch_str}/test.dat already exists -> nothing to be done"
      else
         echo "read_binary_station_beam_test_order1_2pol $dat_file  -f ${ch_str}/test -p 0 -C ${n_channels} -c ${ch} -s ${start_byte} -Z \"$conversion_options\""
         read_binary_station_beam_test_order1_2pol $dat_file  -f ${ch_str}/test -p 0 -C ${n_channels} -c ${ch} -s ${start_byte} -Z "$conversion_options"
      fi
      
      if [[ $process_channel -gt 0 ]]; then
            cd ${ch_str}
            pwd
            if [[ -s ../header_template.txt ]]; then
               # awk -v channel=${ch_full} '{if($1=="FREQ"){printf("FREQ %.4f\n",channel*(400.00/512.00));}else{print $0;}}' ../header_template.txt | head --bytes=4096 > header.txt
               echo "python $path_dada_converter ../header_template.txt header.txt --channel=${ch_full}"
               python $path_dada_converter ../header_template.txt header.txt --channel=${ch_full}      
            
               header_size=`du -sb header.txt | awk '{print $1;}'`
               echo "INFO : header size = $header_size bytes"
               if [[ $header_size != 4096 ]]; then
                  echo "ERROR : .dada header size should be 4096 bytes and is $header_size bytes -> exiting now"
                  exit
               fi
            
               echo "cat header.txt test.dat > test.dada"
               cat header.txt test.dat > test.dada
                        
               if [[ $do_psr_processing -gt 0 ]]; then
                  echo "INFO : performing pulsar processing"
                  
                  if [[ -s ../${object}.eph ]]; then                              
                     echo "cp ../${object}.eph ."
                     cp ../${object}.eph .
                  else
                     echo "psrcat -e ${object} > ${object}.eph"
                     psrcat -e ${object} > ${object}.eph                              
                  fi
              
                  if [[ -s ${object}.eph ]]; then
                     echo "OK : ephemeris file ${object}.eph exists -> processing continues"
                  else
                     echo "ERROR : ephemeris file ../${object}.eph not found and could not be generated -> exiting"
                     exit
                  fi
            
                  size_mb=`du -sm test.dat | awk '{print $1;}'`
            
                  echo "dspsr -E ${object}.eph -b 64 -U ${size_mb} -F 256:D -F 256:D test.dada"
                  dspsr -E ${object}.eph -b 64 -U ${size_mb} -F 256:D -F 256:D test.dada
            
                  echo "rm -f test.dat"
                  rm -f test.dat
               else
                  pwd
                  echo "ERROR : .dada header template ../header_template.txt not found -> cannot process single channels"
               fi
            else
               echo "WARNING : pulsar processing is not required"
            fi
            cd ..
      else 
         echo "WARNING : processing of channels is not required"
      fi
      
      if [[ $test_pattern_generator -gt 0 ]]; then
         # test_pattern_generated_data channel_0_16_1644918011.076890.dat -C 16  > test.out 2>&1
         cd ${ch_str}
         start_value=$((${ch}*4+1))         
         echo "test_pattern_generated_data test.dat -C 1 -S $start_value > test.out 2>&1"
         test_pattern_generated_data test.dat -C 1 -S $start_value > test.out 2>&1
         cd ..
      else
         echo "INFO : testing of pattern generated in firmware is not requested"
      fi
   else
      echo "DEBUG : channel ${ch_str} already processed"
   fi
   
   ch=$(($ch+1))
done


quick_fine_chan_subdir.sh test.dat
merge_channels.sh normal -n
