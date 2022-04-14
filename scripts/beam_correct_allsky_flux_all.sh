#!/bin/bash

for logfile in `ls *-xx.txt  ls *-yy.txt`
do
   echo "beam_correct_allsky_flux.sh ${logfile}"
   beam_correct_allsky_flux.sh ${logfile}
done
