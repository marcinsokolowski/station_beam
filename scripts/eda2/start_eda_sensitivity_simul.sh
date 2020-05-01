#!/bin/bash

channel=204
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

path=`pwd`
dir_path=`dirname $path`
utc=`basename $dir_path`

ux_start=`date2date -ut2ux=$utc | awk '{print $3;}'`
mwa_coarse_channel=`echo $channel | awk '{mwa_cc=$1*((400.00/512.00)/1.28);printf("%d\n",mwa_cc);}'` # this is MWA 1.28-MHz coarse channel ( 107*1.28 = 136.96 MHz )

export PATH=~/github/station_beam/:/opt/anaconda2/bin:$PATH

echo "nohup calc_eda_sensitvity_all.sh $ux_start 86400 - \"--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2\" $mwa_coarse_channel - 300 400 400 > out 2>&1 &"
nohup calc_eda_sensitvity_all.sh $ux_start 86400 - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" $mwa_coarse_channel - 300 400 400 > out 2>&1 &


# ~/github/station_beam/tools/merge_sensitivity_loop.sh