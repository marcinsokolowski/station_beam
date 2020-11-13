#!/bin/bash

# generate_beam_on_sun_file.sh

station=eda2
if [[ -n "$1" && "$1" != "-" ]]; then
   station=$1
fi

beams_name=EDA
if [[ $station != "eda2" ]]; then
   beams_name=AAVS2
fi

local_dir=/media/msok/0754e982-0adb-4e33-80cc-f81dda1580c8/${station}/data/real_time_calibration/
cd ${local_dir}
pwd

dtm=2018_01_01-00:00
if [[ -n "$2" && "$2" != "-" ]]; then
   dtm=$2
else
   path=`ssh ${station} "ls -dtr /data/real_time_calibration/????_??_??-??:??" | tail -1`
   dtm=`basename $path`
fi

echo "rsync --exclude 'cal.out' --exclude 'NoSunBeamCorr' -avP ${station}:/data/real_time_calibration/${dtm} ."
rsync --exclude 'cal.out' --exclude 'NoSunBeamCorr' -avP ${station}:/data/real_time_calibration/${dtm} .

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

