#!/bin/bash

beam_dir=~/aavs-calibration/
update_config=0
if [[ -n "$1" && "$1" != "-" ]]; then
   beam_dir=$1
   update_config=1
fi

beam_url_eda="https://www.dropbox.com/s/q4v9kuf7x1avmvw/EDA.tar.gz?dl=0"
if [[ -n "$2" && "$2" != "-" ]]; then
   beam_url_eda=$2
fi

beam_url_skala2="https://www.dropbox.com/s/vhoedc5geulkeus/SKALA2.tar.gz?dl=0"
if [[ -n "$3" && "$3" != "-" ]]; then
   beam_url_skala2=$3
fi

beam_url_skala4="https://www.dropbox.com/s/n5oz2kc1tejxo2s/SKALA4.tar.gz?dl=0"
if [[ -n "$4" && "$4" != "-" ]]; then
   beam_url_skala4=$4
fi

echo "############################################################################"
echo "PARAMETERS:"
echo "############################################################################"
echo "beam_dir        = $beam_dir"
echo "beam_url_eda    = $beam_url_eda"
echo "beam_url_skala2 = $beam_url_skala2"
echo "beam_url_skala4 = $beam_url_skala4"
echo "############################################################################"


mkdir -p /tmp/station_beam/BeamModels/
cd /tmp/station_beam/BeamModels/

echo "rm -f EDA.tar.gz SKALA.tar.gz SKALA4.tar.gz"
rm -f EDA.tar.gz SKALA.tar.gz SKALA4.tar.gz 

echo "wget ${beam_url_eda}"
wget ${beam_url_eda}

echo "wget ${beam_url_skala2}"
wget ${beam_url_skala2}

echo "wget ${beam_url_skala4}"
wget ${beam_url_skala4}

echo "tar zxvf EDA.tar.gz"
tar zxvf EDA.tar.gz
ln -s EDA EDA2

echo "tar zxvf SKALA2.tar.gz"
tar zxvf SKALA2.tar.gz

echo "tar zxvf SKALA4.tar.gz"
tar zxvf SKALA4.tar.gz

echo "rm -f EDA.tar.gz SKALA.tar.gz SKALA4.tar.gz"
rm -f EDA.tar.gz SKALA.tar.gz SKALA4.tar.gz 

# echo "tar zxvf BeamModels.tar.gz"
# tar zxvf BeamModels.tar.gz

cd ../
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




