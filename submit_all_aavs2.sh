#!/bin/bash

za=0


# ./pawsey_start_calc_eda_sensitvity_timestep_eda2.sh 1578441600
ux=1578441600
end_ux=$(($ux+86400))
ux_step=1800

while [[ $ux -le $end_ux ]];
do
   echo "./pawsey_start_calc_eda_sensitvity_timestep_aavs2.sh ${ux}"
   ./pawsey_start_calc_eda_sensitvity_timestep_aavs2.sh ${ux}
   
   echo "sleep 5"
   sleep 5
   
   ux=$(($ux+$ux_step))
done
