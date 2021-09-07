#!/bin/bash

sql_dir=~/github/station_beam/sql/
if [[ -n "$1" && "$1" != "-" ]]; then
   sql_dir=$1
fi

url_eda="https://www.dropbox.com/s/r2q2gwci18opg3l/ska_station_sensitivity_EDA2.db.gz?dl=0"
if [[ -n "$2" && "$2" != "-" ]]; then
   url_eda=$2
fi

url_aavs2="https://www.dropbox.com/s/vxm407brvbo0e6g/ska_station_sensitivity_AAVS2.db.gz?dl=0"
if [[ -n "$3" && "$3" != "-" ]]; then
   url_aavs2=$3
fi


echo "mkdir -p ${sql_dir}"
mkdir -p ${sql_dir}

cd ${sql_dir}/
echo "wget ${url_eda} -O ska_station_sensitivity_EDA2.db.gz"
wget ${url_eda} -O ska_station_sensitivity_EDA2.db.gz

echo "gzip -d ska_station_sensitivity_EDA2.db.gz"
gzip -d ska_station_sensitivity_EDA2.db.gz


echo "wget ${url_aavs2} -O ska_station_sensitivity_AAVS2.db.gz"
wget ${url_aavs2} -O ska_station_sensitivity_AAVS2.db.gz

echo "gzip -d ska_station_sensitivity_AAVS2.db.gz"
gzip -d ska_station_sensitivity_AAVS2.db.gz

