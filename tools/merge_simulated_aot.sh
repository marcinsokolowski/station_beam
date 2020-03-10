#!/bin/bash

# cat ??????????_az0_za0_eda_sensitivity_XX.txt | grep -v "#"

pol=XX
if [[ -n "$1" && "$1" != "-" ]]; then
   pol=$1
fi

for file in `ls ??????????_az0_za0_eda_sensitivity_${pol}.txt`
do
   gps=`echo $file | cut -b 1-10`
   ux=`gps2ux! $gps`
   aot=`cat $file | grep -v "#" | tail -1 | awk '{print $2;}'`

   if [[ -n $aot ]]; then   
      echo $ux" "$aot
   fi
done
