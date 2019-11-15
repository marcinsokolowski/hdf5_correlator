#!/bin/bash

options=""
if [[ -n "$1" && "$1" != "-" ]]; then
   options="$1"
fi

outfile=beamformed2.txt
if [[ -n "$2" && "$2" != "-" ]]; then
   outfile=$2
fi

do_phase_calib=1
if [[ -n "$3" && "$3" != "-" ]]; then
   do_phase_calib=$3
fi

echo "rm -f ${outfile}"
rm -f ${outfile}

#!/bin/bash
# was : for hdf5file in `ls -tr channel_cont_*.hdf5`
for hdf5file in `ls channel_cont_*.hdf5`
do
   out_txtfile=${hdf5file%%hdf5}txt
   # echo "hdf5_correlator ${hdf5file} -B -o beamformed2.txt -f -x 1000"
   # hdf5_correlator ${hdf5file} -B -o beamformed2.txt -f -x 1000  

   echo "hdf5_correlator ${hdf5file} -B -o ${outfile} -f -c ${do_phase_calib} ${options}"
   hdf5_correlator ${hdf5file} -B -o ${outfile} -f -c ${do_phase_calib} ${options}
done


out_power_file=${outfile%%.txt}_power_vs_time.txt
awk '{print $6" "$5;}' ${outfile} > ${out_power_file}
