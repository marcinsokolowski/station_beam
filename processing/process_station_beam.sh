#!/bin/bash

station_file=`ls *.hdf5  |tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   station_file=$1
fi
channel=4

beam_scripts_path=~/github/station_beam/processing/

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 > x.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 > x.out 2>&1

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 > y.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 > y.out 2>&1

# root -l "plot_vs_time.C(\"power_vs_time_ch${channel}.txt\",-1e6,1e6,NULL,0,1000)"

