#!/bin/bash

# cd /data/msok/sensitivity/db/EDA2/za0-30deg
# ls -trl *.txt | awk '{print substr($9,1,10);}' | sort -u -n

station=EDA2
if [[ -n "$1" && "$1" != "-" ]]; then
   station=$1
fi

cd /home/msok/github/station_beam/data/${station}/nimbus5

while [ 1 ];
do
   rsync -avP --exclude='*.fits' --exclude='*out' nimbus5:/data/msok/sensitivity/db/${station}/za* .

   cd /home/msok/github/station_beam/sql
   ./load_all.sh /home/msok/github/station_beam/data/${station}/nimbus5/za0-30deg - ${station}
   ./load_all.sh /home/msok/github/station_beam/data/${station}/nimbus5/za30-60deg - ${station}
   ./load_all.sh /home/msok/github/station_beam/data/${station}/nimbus5/za60-90deg - ${station}
   cd -

   sleep 300
done
   

