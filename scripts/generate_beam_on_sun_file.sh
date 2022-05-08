#!/bin/bash

station_name=EDA
if [[ -n "$1" && "$1" != "-" ]]; then
   station_name=$1
fi

ls -tr *.hdf5 > hdf5_list

path=`which fits_beam.py`

echo "python $path --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt --station=${station_name}"
python $path --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt --station=${station_name}

