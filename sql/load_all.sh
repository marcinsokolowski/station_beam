#!/bin/bash

data_dir=/home/msok/github/station_beam/data/EDA2/nimbus5/za0-30deg
if [[ -n "$1" && "$1" != "-" ]]; then
   data_dir=$1
fi

do_convert=1
if [[ -n "$2" && "$2" != "-" ]]; then
   do_convert=$2  
fi

station_name=EDA2
if [[ -n "$3" && "$3" != "-" ]]; then
   station_name=$3
fi

station_beam_dir=~/github/station_beam/
if [[ -n "$4" && "$4" != "-" ]]; then
   station_beam_dir=$4
fi

do_load=0
if [[ -n "$5" && "$5" != "-" ]]; then
   do_load=$5
fi



cd ${data_dir}/
rm -f all.sql
for txtfile in `find . -name "1*.txt"`
do
   sqlfile=${txtfile%%txt}sql
 
   echo "${station_beam_dir}/sql/txt2sql.sh $txtfile $sqlfile"
   ${station_beam_dir}/sql/txt2sql.sh $txtfile $sqlfile
   
   echo "cat $sqlfile >> all.sql"
   cat $sqlfile >> all.sql
done


if [[ $do_load -gt 0 ]]; then
   date
   echo "sqlite3 ${station_beam_dir}/sql/ska_station_sensitivity_${station_name}.db <  all.sql"
   sqlite3 ${station_beam_dir}/sql/ska_station_sensitivity_${station_name}.db <  all.sql 
   date
else
   date
   echo "WARNING : loading to the database is not required"
fi   
   
# for txtfile in `ls 1*.txt`
#do
#   sqlfile=${txtfile%%txt}sql
#      
#   echo
#   date
#   echo "sqlite3 ${station_beam_dir}/sql/ska_station_sensitivity_${station_name}.db < ${sqlfile}"
#   sqlite3 ${station_beam_dir}/sql/ska_station_sensitivity_${station_name}.db < ${sqlfile}   
#done
#cd -
