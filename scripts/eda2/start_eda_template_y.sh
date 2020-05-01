#!/bin/bash

# sleep 86400
# sleep 86400
# sleep 86400
# sleep 86400


export PATH="/opt/anaconda2/bin:$PATH"

start_date=20200409_170100
n_days=1
pointing_az_deg=0.00
pointing_za_deg=0.00
time_step=60
parallel=1
freq=160 # ch=204 -> 159.375

nohup simulate_total_power_vs_time_date_range.sh $start_date $n_days $time_step /raid/data/eda/eda2/simulations /raid/data/eda/eda2/simulations/AntennaPatterns/EDA2_256/EDA2_256_elem_station_${freq}MHz_Ypol-DONT-DO.txt /home/bighorns/bighorns/software/analysis/scripts/shell/eda/eda2/config/antenna_locations_eda2.txt $parallel ${freq} - "YY" ${pointing_az_deg} ${pointing_za_deg} > out 2>&1 &
echo "nohup simulate_total_power_vs_time_date_range.sh $start_date $n_days $time_step /raid/data/eda/eda2/simulations /raid/data/eda/eda2/simulations/AntennaPatterns/EDA2_256/EDA2_256_elem_station_${freq}MHz_Ypol-DONT-DO.txt /home/bighorns/bighorns/software/analysis/scripts/shell/eda/eda2/config/antenna_locations_eda2.txt $parallel ${freq} - "YY" ${pointing_az_deg} ${pointing_za_deg} > out 2>&1 &"

