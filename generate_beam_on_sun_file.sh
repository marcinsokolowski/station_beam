#!/bin/bash

ls *.hdf5 > hdf5_list

echo "python ~/github/station_beam/fits_beam.py --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt"
python ~/github/station_beam/fits_beam.py --infile_hdf5list=hdf5_list --outfile_beam_on_sun=beam_on_sun.txt

