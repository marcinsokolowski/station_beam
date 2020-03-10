#!/bin/bash

sleep_time=10
if [[ -n "$1" && "$1" != "-" ]]; then
   sleep_time=$1
fi

while [ 1 ];
do
   merge_simulated_aot.sh XX > aot_vs_ux_X.txt
   merge_simulated_aot.sh YY > aot_vs_ux_Y.txt

   aot2sefd.sh aot_vs_ux_X.txt > sefd_vs_ux_X.txt
   aot2sefd.sh aot_vs_ux_Y.txt > sefd_vs_ux_Y.txt

   echo "sleep $sleep_time"
   sleep $sleep_time
done
