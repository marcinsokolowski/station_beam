#!/bin/bash

station_beam_dir=/home/msok/github/station_beam/

data_dir=/home/msok/github/station_beam/data/EDA2/nimbus5/za0-30deg
if [[ -n "$1" && "$1" != "-" ]]; then
   data_dir=$1
fi


${station_beam_dir}/sql/txt2sql.sh

for sqlfile in `ls ${data_dir}/*.sql`
do
   echo
   date
   echo "sqlite3 ska_station_sensitivity.db < ${sqlfile}"
   sqlite3 ska_station_sensitivity.db < ${sqlfile}   
done


