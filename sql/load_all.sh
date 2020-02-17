#!/bin/bash

station_beam_dir=/home/msok/github/station_beam/

data_dir=/home/msok/github/station_beam/data/EDA2/nimbus5/za0-30deg
if [[ -n "$1" && "$1" != "-" ]]; then
   data_dir=$1
fi

do_convert=1
if [[ -n "$2" && "$2" != "-" ]]; then
   do_convert=$2  
fi

if [[ $do_convert -gt 0 ]]; then
   cd ${data_dir}/
   echo "${station_beam_dir}/sql/txt2sql.sh"
   ${station_beam_dir}/sql/txt2sql.sh
   cd -
else 
   echo "WARNING : conversion step not required (set 2nd parameter to 1 to enforce)"
fi   

for sqlfile in `ls ${data_dir}/*.sql`
do
   echo
   date
   echo "sqlite3 ska_station_sensitivity.db < ${sqlfile}"
   sqlite3 ska_station_sensitivity.db < ${sqlfile}   
done



