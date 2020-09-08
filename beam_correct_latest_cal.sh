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

path=`ssh ${station} "ls -dtr /data/real_time_calibration/????_??_??-??:??" | tail -1`
dtm=`basename $path`

echo "rsync -avP ${station}:/data/real_time_calibration/${dtm} ."
rsync -avP ${station}:/data/real_time_calibration/${dtm} .

if [[ -d ${dtm} ]]; then
   cd ${dtm}
   
   echo "generate_beam_on_sun_file.sh ${beams_name} > beam_on_sun.out 2>&1"
   generate_beam_on_sun_file.sh ${beams_name} > beam_on_sun.out 2>&1
   
   # copy to server :
   ssh ${station} "mkdir -p /tmp/msok"
   
   echo "scp beam_on_sun.txt ${station}:/tmp/msok"
   scp beam_on_sun.txt ${station}:/tmp/msok
   
   cd -
else
   echo "ERROR : directory ${dtm} does not exist"
fi

