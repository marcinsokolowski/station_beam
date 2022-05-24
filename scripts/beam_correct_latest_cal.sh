#!/bin/bash

# generate_beam_on_sun_file.sh

# add paths:
export PATH=/home/aavs/Software/station_beam/scripts/:/home/aavs/Software/station_beam/python/:$PATH

station=eda2
if [[ -n "$1" && "$1" != "-" ]]; then
   station=$1
fi

beams_name=EDA
if [[ $station != "eda2" ]]; then
   beams_name=AAVS2
fi

local_dir=`pwd`
if [[ -n "$3" && "$3" != "-" ]]; then
   local_dir=$3
fi

cd ${local_dir}
pwd

dtm=2018_01_01-00:00
if [[ -n "$2" && "$2" != "-" ]]; then
   dtm=$2
else
   path=`ssh ${station} "ls -dtr /data/real_time_calibration/????_??_??-??:??" | tail -1`
   dtm=`basename $path`
fi

do_copy=0
if [[ -n "$4" && "$4" != "-" ]]; then
   do_copy=$4
fi

echo "####################################################################"
echo "PARAMETERS:"
echo "####################################################################"
echo "station = $station"
echo "beam_name = $beam_name"
echo "local_dir = $local_dir"
echo "do_copy = $do_copy"
echo "####################################################################"
date


if [[ $do_copy -gt 0 ]]; then
   echo "rsync --exclude 'cal.out' --exclude 'NoSunBeamCorr' -avP ${station}:/data/real_time_calibration/${dtm} ."
   rsync --exclude 'cal.out' --exclude 'NoSunBeamCorr' -avP ${station}:/data/real_time_calibration/${dtm} .
else
   echo "WARNING : data copying is not required"
fi   

if [[ -d ${local_dir}/${dtm} ]]; then
   cd ${local_dir}/${dtm}
   
   echo "gzip -df *.hdf5.gz"
   gzip -df *.hdf5.gz
   
   echo "generate_beam_on_sun_file.sh ${beams_name} > beam_on_sun.out 2>&1"
   generate_beam_on_sun_file.sh ${beams_name} > beam_on_sun.out 2>&1
   
   # copy to server :
   ssh ${station} "mkdir -p /tmp/msok"
   
   echo "scp beam_on_sun.txt ${station}:/tmp/msok"
   scp beam_on_sun.txt ${station}:/tmp/msok
   
   cd -
else
   echo "ERROR : directory ${local_dir}/${dtm} does not exist"
fi

