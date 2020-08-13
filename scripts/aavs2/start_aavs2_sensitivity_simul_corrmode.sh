#!/bin/bash

channel=204
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

interval=86400
if [[ -n "$2" && "$2" != "-" ]]; then
   interval=$2
fi

options="" # --correlated_mode
if [[ -n "$3" && "$3" != "-" ]]; then
   options=$3
fi


path=`pwd`
dir_path=`dirname $path`

local_time=`basename $dir_path`
ux_start=`date2date -local2ux=$local_time | awk '{print $1;}'`

# mwa_coarse_channel=`echo $channel | awk '{mwa_cc=$1*((400.00/512.00)/1.28);printf("%d\n",mwa_cc);}'` # this is MWA 1.28-MHz coarse channel ( 107*1.28 = 136.96 MHz )
mwa_coarse_channel=`echo $channel | awk '
function round(x)
{
   ival = int(x)    # integer part, int() truncates

   # see if fractional part
   if (ival == x)   # no fraction
      return ival   # ensure no decimals

   if (x < 0) {
      aval = -x     # absolute value
      ival = int(aval)
      fraction = aval - ival
      if (fraction >= .5)
         return int(x) - 1   # -2.5 --> -3
      else
         return int(x)       # -2.3 --> -2
   } else {
      fraction = x - ival
      if (fraction >= .5)
         return ival + 1
      else
         return ival
   }
}
{
	mwa_cc=$1*((400.00/512.00)/1.28);print round(mwa_cc);
}'`


export PATH=~/github/station_beam/:/opt/anaconda2/bin:$PATH

# echo "nohup calc_eda_sensitvity_all.sh $ux_start 86400 - \"--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2\" $mwa_coarse_channel - 300 400 400 > out 2>&1 &"
# nohup calc_eda_sensitvity_all.sh $ux_start 86400 - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" $mwa_coarse_channel - 300 400 400 > out 2>&1 &

if [[ ! -s antenna_locations_aavs2.txt ]]; then
   echo "cp ~/aavs-calibration/config/aavs2/antenna_locations.txt antenna_locations_aavs2.txt"
   cp ~/aavs-calibration/config/aavs2/antenna_locations.txt antenna_locations_aavs2.txt
else
   echo "File antenna_locations_aavs2.txt already exist -> copy skipped"
fi

# --correlated_mode
echo "nohup calc_eda_sensitvity_all.sh $ux_start ${interval} - \"--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt --correlated_mode ${options}\" $mwa_coarse_channel - 300 400 400 > out 2>&1 &"
nohup calc_eda_sensitvity_all.sh $ux_start ${interval} - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt --correlated_mode ${options}" $mwa_coarse_channel - 300 400 400 > out 2>&1 &

# ~/github/station_beam/tools/merge_sensitivity_loop.sh