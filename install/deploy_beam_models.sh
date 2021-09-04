#!/bin/bash

beam_dir=~/aavs-calibration/
update_config=0
if [[ -n "$1" && "$1" != "-" ]]; then
   beam_dir=$1
   update_config=1
fi

beam_url=
if [[ -n "$2" && "$2" != "-" ]]; then
   beam_url=$2
fi

mkdir -p /tmp/station_beam
cd /tmp/station_beam
echo "wget ${beam_url}"
wget ${beam_url}

echo "tar zxvf BeamModels.tar.gz"
tar zxvf BeamModels.tar.gz

echo "mv BeamModels ${beam_dir}/"
mv BeamModels ${beam_dir}/

if [[ $update_config -gt 0 ]]; then
   # update config file :
   awk -v beam_dir=${beam_dir} -F "=" '{if($1=="beam_model_path"){print $1"="beam_dir;}else{print $0;}}' ~/github/station_beam/station_beam_config.py > /tmp/station_beam/station_beam_config.py

   echo "cp ~/github/station_beam/station_beam_config.py ~/github/station_beam/station_beam_config.py.backup"
   cp ~/github/station_beam/station_beam_config.py ~/github/station_beam/station_beam_config.py.backup

   echo "cp /tmp/station_beam/station_beam_config.py ~/github/station_beam/station_beam_config.py"
   cp /tmp/station_beam/station_beam_config.py ~/github/station_beam/station_beam_config.py
else
   echo "INFO : update of config file ~/github/station_beam/station_beam_config.py is not required (no change to default path to beam directory = $beam_dir)"
fi   




