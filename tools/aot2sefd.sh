#!/bin/bash


file=aot_vs_ux_X.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   file=$1
fi

awk '{
   ux=$1;
   aot=$2;
   sefd=(2.00*1380.00)/aot;
   printf("%.1f %.8f %.8f\n",ux,sefd,aot);
}' ${file}
