#!/bin/bash

# ./pawsey_start_calc_eda_sensitvity_timestep_eda2.sh 1578441600
ux=1578441600
if [[ -n "$1" && "$1" != "-" ]]; then
   ux=$1
fi

interval=86400
if [[ -n "$2" && "$2" != "-" ]]; then
   interval=$2
fi


end_ux=$(($ux+$interval))
ux_step=1800

while [[ $ux -le $end_ux ]];
do
   echo "./pawsey_start_calc_eda_sensitvity_timestep_eda2.sh ${ux}"
   ./pawsey_start_calc_eda_sensitvity_timestep_eda2.sh ${ux}
   
   echo "sleep 5"
   sleep 5
   
   ux=$(($ux+$ux_step))
done
