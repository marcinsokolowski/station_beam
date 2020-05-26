#!/bin/bash

channel=204
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi


path=`pwd`
dir_path=`dirname $path`

local_time=`basename $dir_path`
ux_start=`date2date -local2ux=$local_time | awk '{print $1;}'`

mwa_coarse_channel=`echo $channel | awk '{mwa_cc=$1*((400.00/512.00)/1.28);printf("%d\n",mwa_cc);}'`

export PATH=~/github/station_beam/:/opt/anaconda2/bin:$PATH

# echo "nohup calc_eda_sensitvity_all.sh $ux_start 86400 - \"--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2\" $mwa_coarse_channel - 300 400 400 > out 2>&1 &"
# nohup calc_eda_sensitvity_all.sh $ux_start 86400 - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" $mwa_coarse_channel - 300 400 400 > out 2>&1 &

echo "cp ~/aavs-calibration/config/aavs2/antenna_locations.txt antenna_locations_aavs2.txt"
cp ~/aavs-calibration/config/aavs2/antenna_locations.txt antenna_locations_aavs2.txt

echo "nohup calc_eda_sensitvity_all.sh $ux_start 86400 - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt" $mwa_coarse_channel - 300 400 400 > out 2>&1 &"
nohup calc_eda_sensitvity_all.sh $ux_start 86400 - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt" $mwa_coarse_channel - 300 400 400 > out 2>&1 &

# ~/github/station_beam/tools/merge_sensitivity_loop.sh