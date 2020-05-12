#!/bin/bash

path=`which hdf5_to_dada_converter.py`

for datfile in `ls *.dat`
do
   echo "python $path ${datfile} --dat2dada"
   python $path ${datfile} --dat2dada
done
