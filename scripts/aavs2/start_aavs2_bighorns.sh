#!/bin/bash

export PATH="/opt/anaconda_new/bin:$PATH"
export PATH=~/github/station_beam/:$PATH


# start nimbus5 t0=1578441600 + 3 hours = 1578441600 + 3*3600 = 1578452400 

start_time=1578452400
interval=10800 # 3 hours 


za_start=0 # 30
za_end=25
nohup calc_eda_sensitvity_all.sh $start_time $interval - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt" - - 1800 - - - ${za_start} ${za_end} > za0-25deg.out 2>&1 &

sleep 10
za_start=30 # 30
za_end=55
nohup calc_eda_sensitvity_all.sh $start_time $interval - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt" - - 1800 - - - ${za_start} ${za_end} > za30-55deg.out 2>&1 &

sleep 10
za_start=60 # 30
za_end=90
nohup calc_eda_sensitvity_all.sh $start_time $interval - "--use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_aavs2.txt" - - 1800 - - - ${za_start} ${za_end} > za60-90deg.out 2>&1 &
