#!/bin/bash

do_correlation=0
bin2lfiles=0
convert2casa=1
dumpbinfile=1
ra_hrs=
dec_degs=
# lfile_converter_path="./"
# (1.13246208*1000000.00)/(32*1.08) = 32768
n_avg=32700 # ~1 second and 2894 was ok for 0.1 second 
n_chan=32
inttime=1.13246208
corr_path="/home/rwayth/bin/corr_gpu_complex2"
aavs_calibration_path=~/aavs-calibration/ # WAS : "/home/rwayth/aavs-calibration/" # Lfile2uvfits_eda.sh -i 0.1 -n 11 -R 8.19680715619 -D 19.9892573884 -N 512 -C 32 20190724_041444_eda2_ch32_ant256_midday_avg2894
force=0
auto_sun=1 # automatically set phase center to Sun
timestamp_from_hdf5file=1
merged_dir="merged/"
hdf5_file_list="hdf5_file_list.tmp"
n_integrations_per_uvfits=1
n_samples_per_hdf5_file=-1
freq_channel=204
freq_channel_param=-1
do_merge=1
remove_single_hdf5_files=0
antenna_locations_path="antenna_locations.txt"
instr_path="instr_config.txt"
dtm_ux=0
hdf5_template="channel_cont_0_????????_*_*.hdf5"
channelised_data=1
station_name="eda2"

if [[ -s /opt/aavs/config/station.yml ]]; then
   station_name=`awk -v station_section=0 '{if(index($1,":")>0 && NF==1){if(index($1,"station")>0 ){station_section=1;}else{station_section=0;}}if(station_section>0){if($1=="name:"){station_name=$2;gsub("\"","",station_name);print tolower(station_name);}}}' /opt/aavs/config/station.yml`   
   echo "Station config file (or symbolik link) exists -> getting station_name = $station_name"
else
   echo "WARNING : /opt/aavs/config/station.yml file or symbolic link does not exist will use default station_name = $station_name or value passed in parameter -s"   
fi

function print_usage {
  echo "Script merges all hdf5 files from 16 tiles to a single one, converts them into bin files and depending on options runs correlation and lfiles to uvfits conversion"
  echo "Usage: "
  echo "hdf5_to_uvfits_all.sh [options]"
  echo "    -c Enables correlation Default: $do_correlation"
  echo "    -C CORRELATOR_PROGRAM_PATH , default $corr_path"
  echo "    -l Converts l-files to uvfits files Default: $bin2lfiles"
  echo "    -R ra_hours   Default: use zenith"
  echo "    -D dec_degs   Default: use zenith"
  echo "    -F forces overwritting of bin file"
  echo "    -i INTTIME ( calculate as number of samples collected typically 131072 x 1.08 usec  = 0.14155776 sec or when I did 1048576 samples it is 1.13246208 sec [default $inttime]"
  echo "    -n n_chunks (same as -n n_chunks in Randall's Lfile2uvfits_eda.sh) to have 1 integration per uvfits for n samples 1048576 use 32700 and for 131072 samples use 4070 (n_samples/32 fft samples) [default $n_avg]"
  echo "    -z : turn on automating setting of phase center to zenith [default auto_sun=$auto_sun]"
  echo "    -t : get timestamps from hdf5 sample_timestamps table [default $timestamp_from_hdf5file]"
  echo "    -d merged_directory : where to save or expect to find merged .hdf5 files and .bin files [default $merged_dir]"
  echo "    -L hdf5_file_list : file with list of hdf5 files to convert [default $hdf5_file_list]" 
  echo "    -I n_integrations_per_uvfits : number of integrations per uvfits file [default $n_integrations_per_uvfits]"
  echo "    -H disables dumping of bin files, just merges hdf5 files [default dumpbinfile = $dumpbinfile]"
  echo "    -f frequency channel [default $freq_channel_param]"
  echo "    -N no merge [default do_merge = $do_merge]"
  echo "    -r remove single (original) hdf5 file after merging [default $remove_single_hdf5_files]"
  echo "    -b Number_of_channes [default $n_chan]"
  echo "    -S ENABLE_CONVERT_TO_CASA : enable conversion to CASA measurements set [default $convert2casa]"
  echo "    -T hdf5_template : template of HDF5 file [default $hdf5_template]"
  echo "    -a N_AVG : correlated data from Alessio's correlator [default disabled - assuming channelised data], parameter is number of averages [default $n_avg]"
  echo "    -s STATION_NAME : name of the station"
  echo 
  echo "INFO : -n or -a should correspond to -i if number of samples in voltage dump is 262144 , each 1.08 usec -> 0.2831 seconds -> 262144/32 = 8192 integrations -> -n 8192 -i 0.2831"
  exit
}


# parse command-line args
# if [ $# -lt 1 ] ; then print_usage ; fi
while getopts "HthFclR:D:i:n:zd:L:I:C:f:Nrb:S:T:a:s:" opt; do
  case $opt in
    a)
        # HDF5 files from Alessio correlator 
        channelised_data=0
        hdf5_template="correlation_burst_*_????????_*_*.hdf5"
        do_merge=0
        inttime=1.98180864 # default value corresponding to 1835008 x 1.08 usec / 1000000 = 1.98180864 seconds 
        n_avg=$OPTARG
        n_chan=1 # for now Alessio correlator is just 1 channel 
        ;;
    h)
        print_usage
        ;;
    b)
        n_chan=$OPTARG
        ;;
    S)
        convert2casa=$OPTARG
        ;;
    F)
        force=1
        ;;
    i)
        inttime=$OPTARG
        ;;
    f)
        freq_channel_param=$OPTARG
        tmp=`echo $freq_channel_param | awk '{printf("%d",$1);}'`
        freq_channel_param=$tmp
        ;;
    c)
        do_correlation=1
        ;;
    C)
        corr_path=$OPTARG
        ;;
    l)
        bin2lfiles=1
        ;;
    n)
        n_avg=$OPTARG
        ;;
    R)
        ra_hrs=$OPTARG
        useradec=1
        ;;
    D)
        dec_degs=$OPTARG
        ;;
    d)
        merged_dir=$OPTARG
        ;;
    I)
        n_integrations_per_uvfits=$OPTARG
        ;;
    L)
        hdf5_file_list=$OPTARG
        ;;
    H)
        dumpbinfile=0
        ;;
    N)
        do_merge=0
        ;;
    s)
        station_name=$OPTARG
        ;;
    z)
        auto_sun=0
        ;;
    t)
        timestamp_from_hdf5file=1
        ;;
    r)
        remove_single_hdf5_files=1
        ;;
    T)
        hdf5_template=$OPTARG
        ;;
    \?)
      echo "Invalid option: -$OPTARG" 1>&2
      print_usage
      ;;
  esac
done
shift $(expr $OPTIND - 1 )


echo "##################################################################"
echo "PARAMTERES :"
echo "##################################################################"
echo "station_name   = $station_name"
echo "hdf5_file_list = $hdf5_file_list"
echo "hdf5_template  = $hdf5_template"
echo "channelised_data = $channelised_data (merge = $do_merge)"
echo "freq_channel   = $freq_channel_param"
echo "do_correlation = $do_correlation" 
echo "bin2lfiles     = $bin2lfiles"
echo "convert2casa   = $convert2casa"
echo "dumpbinfile    = $dumpbinfile"
echo "force          = $force"
echo "ra_hrs         = $ra_hrs"
echo "dec_degs       = $dec_degs"
echo "n_avg          = $n_avg"
echo "inttime        = $inttime"
echo "auto_sun       = $auto_sun"
echo "timestamp_from_hdf5file = $timestamp_from_hdf5file"
echo "merged_dir     = $merged_dir"
echo "n_integrations_per_uvfits = $n_integrations_per_uvfits"
echo "n_output_channels = $n_chan"
echo "corr_path      = $corr_path"
echo "do_merge       = $do_merge"
echo "remove_single_hdf5_files = $remove_single_hdf5_files"
echo "##################################################################"


convert_path=`which merge_n_hdf5_files.py`
hdf5_info_path=`which hdf5_info.py`

mkdir -p ${merged_dir}
if [[ $do_correlation -gt 0 ]]; then
   cd ${merged_dir}
   # /home/rwayth/aavs-calibration/antenna_locations_eda2.txt
   station_name_lower=`echo ${station_name} | awk '{print tolower($1);}'`
   
   echo "cp ${aavs_calibration_path}/config/${station_name_lower}/antenna_locations.txt ."
   cp ${aavs_calibration_path}/config/${station_name_lower}/antenna_locations.txt .
   
   
   echo "cp ${aavs_calibration_path}/config/${station_name_lower}/instr_config.txt ."
   cp ${aavs_calibration_path}/config/${station_name_lower}/instr_config.txt .
   
#   echo "cp ${aavs_calibration_path}/config/${station_name_lower}/instr_config_${station_name_lower}.txt ."
#   cp ${aavs_calibration_path}/config/${station_name_lower}/instr_config_${station_name_lower}.txt .
   
#   echo "ln -s instr_config_${station_name_lower}.txt instr_config.txt"
#   ln -s instr_config_${station_name_lower}.txt instr_config.txt

   echo "cp ${aavs_calibration_path}/config/${station_name_lower}/header.txt ."
   cp ${aavs_calibration_path}/config/${station_name_lower}/header.txt .

   cd -
else
   echo "Correlation is not required -> no need to create symbolic links"
fi

echo "Removing automatic file list hdf5_file_list.tmp ( rm -f hdf5_file_list.tmp )"
rm -f hdf5_file_list.tmp

if [[ -s ${hdf5_file_list} ]]; then
   echo "File list $hdf5_file_list already exists"
else
   echo "Creating file $hdf5_file_list as : ls $hdf5_template > $hdf5_file_list"
   ls $hdf5_template > $hdf5_file_list
   count=`cat $hdf5_file_list | wc -l`
   if [[ $count -le 0 ]]; then
      echo "No files channel_cont_0_????????_*_*.hdf5"
      if [[ $do_merge -gt 0 ]]; then
         echo "Setting do_merge := 0"
         do_merge=0
      fi
      ls channel_cont_????????_*_*.hdf5 > $hdf5_file_list
   else
      echo "Non merged files found -> merging"
   fi  
fi

for hdf5_file_tile0 in `cat $hdf5_file_list`
do
    if [[ $freq_channel_param -ge 0 ]]; then
       # user parameter -f overwrites any other frequency requirements :
       freq_channel=$freq_channel_param
    else
       channel_id=`python ~/aavs-calibration/sensitivity/daq/getch.py $hdf5_file_tile0 | grep "is channel" | awk '{print $5;}'`
       echo "channel_id from file $hdf5_file_tile0 is $channel_id"
       freq_channel=$channel_id
    fi

    echo "---------------------------------------------------- $hdf5_file_tile0 (freq_channel = $freq_channel) ----------------------------------------------------"    
    
    lfile_base="Lfile"
    lfile_base_corr=$lfile_base
    if [[ $do_merge -gt 0 ]]; then
       merged_hdf5_file=`echo ${hdf5_file_tile0} | awk '{a=gsub("channel_cont_0_","channel_cont_",$1);print $1;}'`
    else
       merged_hdf5_file=$hdf5_file_tile0
    fi   
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

       echo "INFO FROM hdf5_info.py output : dtm_ux = $dtm_ux (taken from $hdf5_info_file) -> dtm_local = $dtm_local , n_samples_per_hdf5_file = $n_samples_per_hdf5_file"
    else
       dtm_local=`echo $hdf5_file_tile0 | awk '{l=length($1);end=substr($1,l-20);dt=substr(end,1,8);daysec=substr(end,10,5);h=int(daysec/3600.00);m_dec=(((daysec/3600.00)-h)*60.00);m=int(m_dec);s=int((m_dec-int(m))*60.00);printf("%s %02d:%02d:%02d\n",dt,h,m,s);}'`
       # WORKS OK : date -u -d 'TZ="Australia/Perth" 20190724 10:25:07' +%s
       # dtm_string=`echo "'TZ=\"Australia/Perth\" ${dtm_local}'"`
       dtm_ux=`date -d "$dtm_local" +%s`
       echo "dtm_local = $dtm_local -> $dtm_string"
    fi
    dtm_ut=`date -u -d "1970-01-01 UTC $dtm_ux seconds" +"%Y%m%d_%H%M%S"`
    # 20190724_041444_eda2_ch32_ant256_midday_avg2894.bin
    lfile_base="${dtm_ut}_eda2_ch${n_chan}_ant256_midday_avg${n_avg}"
    merged_bin_file=${lfile_base}.bin
 
    radec_string=""   
    if [[ $auto_sun -gt 0  ]]; then # was with : && $bin2lfiles -gt 0
        # old version using my program , new version uses Randall's program 
        # radec_values=`print_sun $dtm_ux -c  | grep "(RA,DEC)"`
        # ra_degs=`echo $radec_values | awk '{print $4;}'`
        # ra_hrs=`echo $radec_values | awk '{print $4/15.00;}'`
        # dec_degs=`echo $radec_values | awk '{print $6;}'`      
        
        dtm_ux_int=`echo $dtm_ux | awk '{print int($1);}'`
        echo "python ${aavs_calibration_path}/sunpos.py $dtm_ux_int"
        radec_values=`python ${aavs_calibration_path}/sunpos.py $dtm_ux_int`
        ra_hrs=`echo $radec_values | awk '{print $1;}'`
        dec_degs=`echo $radec_values | awk '{print $2;}'`
        radec_string="-R $ra_hrs -D $dec_degs"
        
        echo "Sun position : "
        echo "python ${aavs_calibration_path}/sunpos.py $dtm_ux_int"
        python ${aavs_calibration_path}/sunpos.py $dtm_ux_int
    fi

    # command line options can overwrite automatic value
    if [[ -n "$ra_hrs" && -n "$dec_degs" ]]; then
       radec_string="-R $ra_hrs -D $dec_degs"
    fi


    if [[ $channelised_data -gt 0 ]]; then
       echo "$hdf5_file_tile0 : merged_hdf5_file = $merged_hdf5_file , merged_bin_file = $merged_bin_file , $dtm_local, $dtm_ux, $dtm_ut -> lfile_base = $lfile_base , merged_bin_file = $merged_bin_file , radec_string = $radec_string"
       lfile_base_corr=${merged_bin_file%%.bin}
       if [[ ! -s ${merged_dir}/${merged_bin_file} || $force -gt 0 ]]; then
            if [[ ! -s ${merged_dir}/${merged_hdf5_file} || $force -gt 0 ]]; then
#            echo "print_sun $dtm_ux -c"
#            print_sun $dtm_ux -c
#             echo "Sun position : "
#             echo "${aavs_calibration_path}/sunpos.py $dtm_ux"
#             ${aavs_calibration_path}/sunpos.py $dtm_ux
#    exit



        #    echo "hdf5_to_bin.sh $hdf5_tile0_file"
        #    hdf5_to_bin.sh $hdf5_tile0_file       
               echo "python $convert_path ${hdf5_file_tile0} 16 - --outdir=${merged_dir}"
               python $convert_path ${hdf5_file_tile0} 16 - --outdir=${merged_dir}
            
               if [[ $remove_single_hdf5_files -gt 0 ]]; then
                  single_hdf5_template=`echo ${hdf5_file_tile0} | awk '{a=gsub("channel_cont_0_","channel_cont_*_",$1);print $1;}'`
               
                  echo "rm -f ${single_hdf5_template}"
                  rm -f ${single_hdf5_template}
               else
                  echo "INFO : removing of the original single hdf5 files is not required"
               fi
           else
               echo "INFO : ${merged_dir}/${merged_hdf5_file} already exists -> not need to re-create (if want to use -F option to force)"
           fi
        

           if [[ $dumpbinfile -gt 0 ]]; then
               if [[ ! -s ${merged_dir}/${merged_bin_file} || $force -gt 0 ]]; then
                   #    echo "hdf5dump.sh ${merged_hdf5_file}"
                   #    hdf5dump.sh ${merged_hdf5_file}
                   cd ${merged_dir}/
                   pwd    
                   echo "hdf5_correlator ${merged_hdf5_file} -d -o ${merged_bin_file}"
                   hdf5_correlator ${merged_hdf5_file} -d -o ${merged_bin_file}
                   cd -
               else
                   echo "WARNING : ${merged_dir}/${merged_bin_file} already exists -> not need to convert (if want to use -F option to force)"
               fi
           else
               echo "WARNING : dumpbinfile=$dumpbinfile -> conversion from HDF5 to bin is disabled"
           fi
       else
           echo "WARNING : ${merged_dir}/${merged_bin_file} already exists -> not need to merged/convert (if want to use -F option to force)"
       fi

       cd ${merged_dir}    
       if [[ $do_correlation -gt 0 ]]; then
          if [[ -s ${corr_path} ]]; then
             echo "${corr_path} -c ${n_chan} -n 512 -a ${n_avg} -i ${merged_bin_file} -o ${lfile_base} -w 10"
             ${corr_path} -c ${n_chan} -n 512 -a ${n_avg} -i ${merged_bin_file} -o ${lfile_base} -w 10
          
             echo "rm -f ${merged_bin_file}"
             rm -f ${merged_bin_file}
          else
             echo "ERROR : ${corr_path} does not exist -> cannot correlate files"
          fi
       else 
          echo "WARNING : correlation is not required"
       fi 

       if [[ $bin2lfiles -gt 0 ]]; then
          # echo "${lfile_converter_path}/Lfile2uvfits_eda.sh $lfile_base $radec_string"
          # ${lfile_converter_path}/Lfile2uvfits_eda.sh $lfile_base $radec_string
          # NEW : ~/aavs-calibration/Lfile2uvfits_eda.sh -i 0.1 -n 11 -R 8.19680715619 -D 19.9892573884 -N 512 -C ${n_chan} 20190724_041444_eda2_ch${n_chan}_ant256_midday_avg2894

          if [[ -d ${aavs_calibration_path} ]]; then       
             # was -i 1.130112 
             # TEMPORARY unitl Randall commits config files 
             # aavs_calibration_path=~/aavs-calibration/
             echo "${aavs_calibration_path}/Lfile2uvfits_eda.sh -i ${inttime} -n ${n_integrations_per_uvfits} ${radec_string} -N 512 -C ${n_chan} -f ${freq_channel} -F ${lfile_base} -s ${station_name}"
             ${aavs_calibration_path}/Lfile2uvfits_eda.sh -i ${inttime} -n ${n_integrations_per_uvfits} ${radec_string} -N 512 -C ${n_chan} -f ${freq_channel} -F ${lfile_base} -s ${station_name}
          else
             echo "ERROR : ${aavs_calibration_path} does not exist -> cannot convert L-files to uvfits files"
          fi
       else 
          echo "WARNING : conversion from .bin -> lfiles is not required"
       fi 
    else
       # CORRELATED DATA - using single coarse channel xGPU correlator :
       lfile_base_corr=${hdf5_file_tile0%%.hdf5}
       corr_lfile=${lfile_base_corr}.LCCSPC
       auto_lfile=${lfile_base_corr}.LACSPC
    
       if [[ -s $corr_lfile && -s $auto_lfile && $force -le 0 ]]; then
          echo "WARNING : L-files $corr_lfile and $auto_lfile already exist , conversion from .hdf5 to L-files is not required (use option -F to force it)"
       else
          echo "hdf2Lfile.sh ${hdf5_file_tile0} ${n_avg}"
          hdf2Lfile.sh ${hdf5_file_tile0} ${n_avg}
       fi       

       # Modifying integration time accordingly :          
       inttime0=$inttime
       inttime=`echo $inttime $n_avg | awk '{print $1*$2;}'`
       echo "inttime := $inttime ( = $n_avg * $inttime0 )"
       
       if [[ -d ${aavs_calibration_path} ]]; then
           # was -i 1.130112 
           # TEMPORARY unitl Randall commits config files 
           # aavs_calibration_path=~/aavs-calibration/
           echo "${aavs_calibration_path}/Lfile2uvfits_eda.sh -i ${inttime} -n ${n_integrations_per_uvfits} ${radec_string} -N 512 -C ${n_chan} -f ${freq_channel} ${lfile_base_corr} -s ${station_name}"
           ${aavs_calibration_path}/Lfile2uvfits_eda.sh -i ${inttime} -n ${n_integrations_per_uvfits} ${radec_string} -N 512 -C ${n_chan} -f ${freq_channel} ${lfile_base_corr} -s ${station_name}
       else
           echo "ERROR : ${aavs_calibration_path} does not exist -> cannot convert L-files to uvfits files"
       fi
    fi
    
    
    
    if [[ $convert2casa -gt 0 ]]; then
       radec_options=""
       if [[ $auto_sun -gt 0 ]]; then
          ra_deg=`echo $ra_hrs | awk '{print $1*15.00;}'`
          radec_options="-r $ra_deg -d $dec_degs"
       fi
       utc=`date -u -d "1970-01-01 UTC $dtm_ux seconds" +"%Y%m%d_%H%M%S"`
       echo "lfile2casa ${lfile_base_corr} ${lfile_base}.ms -a ${antenna_locations_path} -i ${instr_path} -c ${freq_channel} -u ${utc} ${radec_options} -n ${n_chan} -I ${inttime}"
       lfile2casa ${lfile_base_corr} ${lfile_base}.ms -a ${antenna_locations_path} -i ${instr_path} -c ${freq_channel} -u ${utc} ${radec_options} -n ${n_chan} -I ${inttime}
    else
       echo "WARNING : conversion from L-files to CASA measurements set is not required (use option -S to enable it)"
    fi
    
    cd -
done
