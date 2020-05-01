#!/bin/bash

export PATH="/opt/anaconda2/bin:$PATH"

start_date=20200407_202000
n_days=1
pointing_az_deg=0.00
pointing_za_deg=0.00
time_step=60
parallel=1
freq=160


nohup simulate_total_power_vs_time_date_range.sh $start_date $n_days $time_step /raid/data/aavs/aavs2/simulations /raid/data/eda/eda2/simulations/AntennaPatterns/AAVS2/AAVS2_256_elem_station_${freq}MHz_Xpol.txt.DONT_DO ~/aavs-calibration/config/aavs2/antenna_locations_20191202.txt $parallel ${freq} > out 2>&1 &
echo "nohup simulate_total_power_vs_time_date_range.sh $start_date $n_days $time_step /raid/data/aavs/aavs2/simulations /raid/data/eda/eda2/simulations/AntennaPatterns/AAVS2/AAVS2_256_elem_station_${freq}MHz_Xpol.txt.DONT_DO ~/aavs-calibration/config/aavs2/antenna_locations_20191202.txt $parallel ${freq} > out 2>&1 &"
