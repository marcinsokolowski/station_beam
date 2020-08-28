#!/bin/bash

station_name=EDA
if [[ -n "$1" && "$1" != "-" ]]; then
   station_name=$1
fi

ls -tr *.hdf5 > hdf5_list

echo "python ~/github/station_beam/fits_beam.py --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt --station=${station_name}"
python ~/github/station_beam/fits_beam.py --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt --station=${station_name}

