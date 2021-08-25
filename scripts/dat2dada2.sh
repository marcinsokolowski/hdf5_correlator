#!/bin/bash

object=B0950
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

freq_ch=204
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_ch=$2
fi
freq_ch_val=`echo $freq_ch | awk '{printf("%d",$1);}'`

prefix="*"
if [[ -n "$3" && "$3" != "-" ]]; then
   prefix="$3"
fi

do_dspsr=1
if [[ -n "$4" && "$4" != "-" ]]; then
   do_dspsr=$4
fi

dspsr_options=""
if [[ -n "$5" && "$5" != "-" ]]; then
   dspsr_options=$5
fi

conjugate=1
if [[ -n "$6" && "$6" != "-" ]]; then
   conjugate=$6
fi

force=0
if [[ -n "$7" && "$7" != "-" ]]; then
   force=$7
fi

auto=1
if [[ -n "$8" && "$8" != "-" ]]; then
   auto=$8
fi

phase_bins=64
if [[ -n "$9" && "$9" != "-" ]]; then
   phase_bins=$9
fi


echo "###########################################################"
echo "PARAMETERS:"
echo "###########################################################"
echo "object   = $object"
echo "freq_ch  = $freq_ch"
echo "prefix   = $prefix"
echo "do_dspsr = $do_dspsr"
echo "dspsr_options = $dspsr_options"
echo "conjugate = $conjugate"
echo "force     = $force"
echo "auto      = $auto"
echo "phase_bins = $phase_bins"
echo "###########################################################"



eph_dir=~/github/hdf5_correlator/scripts/config/dspsr/

path=`which hdf5_to_dada_converter.py`

for datfile in `ls ${prefix}.dat`
do
# 2021-08-24 - no longer works after Alessio changed the code to record >1 channel :
#   unixtime=`echo $datfile | cut -b 11-25`
# NEW CODE should handle this new format of dat file name :
   unixtime=`echo $datfile | awk '{l=length($1);ux=substr($1,l-20,17);print ux;}'`
   echo "$datfile -> $unixtime - ok ? (new parsing)"
#   sleep 2

   outfile=${datfile%%.dat}_${object}.dada   
   hdrfile=${datfile%%.dat}_${object}.hdr
  
#   echo "python $path ${datfile} --dat2dada --outfile=${outfile}"
#   python $path ${datfile} --dat2dada --outfile=${outfile}

   if [[ ! -s ${outfile} || $force -gt 0 ]]; then
   
      size_mb=`du -smL ${datfile} | awk '{print $1;}'`
      echo "size_mb = $size_mb"
   
      if [[ $conjugate -gt 0 ]]; then       
         echo "python $path ${datfile} --dat2dada --outfile=${outfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object} --conjugate"
         python $path ${datfile} --dat2dada --outfile=${outfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object} --conjugate
      else
         echo "python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object}"
         python $path ${datfile} --psrdadahdr --outfile=${hdrfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object}
   
         echo "cat ${hdrfile} ${datfile} > ${outfile}"
         cat ${hdrfile} ${datfile} > ${outfile}
#         echo "python $path ${datfile} --dat2dada --outfile=${outfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object}"
#         python $path ${datfile} --dat2dada --outfile=${outfile} --unixtime=${unixtime} --freq_ch=${freq_ch_val} --source=${object}
      fi
   else
      echo "WARNING : dada file ${outfile} already exists -> skipped (enable force=1 in order to re-process)"
   fi            
   
   if [[ $do_dspsr -gt 0 ]]; then
      if [[ ! -s ${object}.eph ]]; then
         if [[ -s ${eph_dir}/${object}.eph ]]; then
            echo "cp ${eph_dir}/${object}.eph ."
            cp ${eph_dir}/${object}.eph .
         else
            echo "psrcat -e ${object} > ${object}.eph"
            psrcat -e ${object} > ${object}.eph
         fi
      fi
   
      if [[ -s ${object}.eph ]]; then
         size_mb=`du -smL ${datfile} | awk '{print $1;}'`
         echo "size_mb = $size_mb"
      
         echo "dspsr -E ${object}.eph -b $phase_bins -U $size_mb -F 256:D ${dspsr_options} ${outfile}"
         dspsr -E ${object}.eph -b $phase_bins -U $size_mb -F 256:D ${dspsr_options} ${outfile}
         
         mkdir -p backup/
         echo "cp *.ar backup/"
         cp *.ar backup/
   
         last_ar=`ls -tr *.ar | tail -1`
   
         echo "psrplot -p flux -D /xs $last_ar"
         psrplot -p flux -D /xs $last_ar
            
         echo "psrplot -p flux -D /png $last_ar"
         psrplot -p flux -D /png $last_ar
            
         pngfile=${last_ar%%ar}png
         echo "mv pgplot.png $pngfile"
         mv pgplot.png $pngfile

         device1=""
         device2=""
         ps_file1=${last_ar%%.ar}_pav1.ps
         ps_file2=${last_ar%%.ar}_pav2.ps
         png_file1=${last_ar%%.ar}_pav1.png
         png_file2=${last_ar%%.ar}_pav2.png
         auto_cmd=""
         if [[ $auto -gt 0 ]]; then
            device1="-g $ps_file1/ps"
            device2="-g $ps_file2/ps"            
            
            auto_cmd="echo | "
         fi
            
         echo "$auto_cmd pav -G -DTp -N1,1 2 $device1 $last_ar"
         $auto_cmd pav -G -DTp -N1,1 2 $device1 $last_ar 
            
         echo "$auto_cmd pav -F -C -d -G -DTp -N1,1 2 $device2 $last_ar"
         pav -F -C -d -G -DTp -N1,1 2 $device2 $last_ar
         
#         if [[ $auto -gt 0 ]]; then
#            echo "convert $ps_file1 $png_file1"
#            convert $ps_file1 $png_file1

#            echo "convert $ps_file2 $png_file2"
#            convert $ps_file2 $png_file2
#         fi
      else
         echo "WARNING : missing file ${object}.eph , cannot find local version neither in ${eph_dir} - please fix it and re-run dspsr"
         echo "dspsr -E ${object}.eph -b 64 -U 600 -F 256:D ${dspsr_options} ${outfile}"
         echo "and : psrplot -p flux -D /xs $last_ar"
      fi
   else
      echo "WARNING : dspsr is not required"
   fi
done

