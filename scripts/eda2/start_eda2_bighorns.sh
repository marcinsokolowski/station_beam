#!/bin/bash

export PATH="/opt/anaconda_new/bin:$PATH"
export PATH=~/github/station_beam/:$PATH


# start nimbus5 t0=1578441600 + 3 hours = 1578441600 + 3*3600 = 1578452400 

start_time=1578452400

za_start=0 # 30
za_end=25
nohup calc_eda_sensitvity_all.sh $start_time - - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" - - 1800 - - - ${za_start} ${za_end} > za0-25deg_out 2>&1 &

za_start=30 # 30
za_end=55
nohup calc_eda_sensitvity_all.sh $start_time - - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" - - 1800 - - - ${za_start} ${za_end} > za30-55deg_out 2>&1 &

za_start=60 # 30
za_end=90
nohup calc_eda_sensitvity_all.sh $start_time - - "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" - - 1800 - - - ${za_start} ${za_end} > za60-90deg_out 2>&1 &

