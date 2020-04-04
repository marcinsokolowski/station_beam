#!/bin/bash

station_file=`ls *.hdf5  |tail -1`
if [[ -n "$1" && "$1" != "-" ]]; then
   station_file=$1
fi
channel=4

tag=`date +%Y%m%d`
if [[ -n "$2" && "$2" != "-" ]]; then
   tag="$2"
fi

last_n_seconds=630720000
if [[ -n "$3" && "$3" != "-" ]]; then
   last_n_seconds=$3
fi

www_dir=aavs1-server:/exports/eda/eda2/station_beam/
if [[ -n "$4" && "$4" != "-" ]]; then
    www_dir=$4
fi

beam_scripts_path=~/github/station_beam/processing/

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} > x.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=0 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} > x.out 2>&1

echo "python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename=\"${tag}_power_vs_time_ch%d_%s.txt\" --last_n_seconds=${last_n_seconds} > y.out 2>&1"
python $beam_scripts_path/hdf5fits_station_beam.py ${station_file} --pol=1 --out_file_basename="${tag}_power_vs_time_ch%d_%s.txt" --last_n_seconds=${last_n_seconds} > y.out 2>&1

ls ${tag}_power_vs_time_ch*.txt > list
first_file=`head --lines=1 list`

mkdir -p images/

root_path=`which root`
if [[ -n $root_path ]]; then
   echo "Plotting in ROOT ..."
   root -l "plotNfiles_vs_time_scaling.C(\"list\",-1,1,-1,1.00)"
else
   echo "WARNING : ROOT cern software not installed, skipping plot in ROOT (no worries !)"
fi

echo "python $beam_scripts_path/plot_power_vs_time.py ${first_file}"
python $beam_scripts_path/plot_power_vs_time.py ${first_file}


png_file=${first_file%%txt}png
# echo "rsync -avP images/${png_file} ${www_dir}/"
# rsync -avP images/${png_file} ${www_dir}/


# root -l "plot_vs_time.C(\"power_vs_time_ch${channel}.txt\",-1e6,1e6,NULL,0,1000)"

