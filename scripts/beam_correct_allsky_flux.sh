#!/bin/bash

input_file=meteor-xx.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   input_file="$1"
fi

freq_mhz=98
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_mhz=$2
fi

pol=`echo $input_file | awk '{l=length($1);pol=substr($1,l-5,2);print toupper(pol);}'`
if [[ -n "$3" && "$3" != "-" ]]; then
   pol=$3
fi

path=`which fits_beam.py`

echo "python $path --time_azh_file=${input_file} --freq_mhz=${freq_mhz} --polarisation=${pol} --projection=\"\""
python $path --time_azh_file=${input_file} --freq_mhz=${freq_mhz} --polarisation=${pol} --projection="" 

echo "python $path --time_azh_file=${input_file} --freq_mhz=${freq_mhz} --polarisation=${pol} --projection=zea"
python $path --time_azh_file=${input_file} --freq_mhz=${freq_mhz} --polarisation=${pol} --projection=zea

