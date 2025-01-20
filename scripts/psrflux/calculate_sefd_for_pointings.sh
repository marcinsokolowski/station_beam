#!/bin/bash

# use anaconda with Python 2.7 :
export PATH="/opt/caastro/ext/anaconda/bin:~/github/station_beam/scripts/psrflux:~/github/station_beam/scripts/$PATH"

infile=data_for_sefd.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   infile="$1"
fi

outfile=${infile%%txt}sefd
if [[ -n "$2" && "$2" != "-" ]]; then
   outfile="$2"
fi

inttime=3600
if [[ -n "$3" && "$3" != "-" ]]; then
   inttime=$3
fi

bw=50
if [[ -n "$4" && "$4" != "-" ]]; then
   bw=$4
fi


if [[ -s ${outfile} ]]; then
   echo "mv ${outfile} ${outfile}.backup"
   mv ${outfile} ${outfile}.backup
fi

options=""

while read line
do
   # psrname RA DEC UNIX Center_Frequency
   # J0034-0534 8.59095042 -5.57683889 1721248417.6704 71.094
   first=`echo $line | awk '{print $1;}'`
   if [[ $first != "#" ]]; then
      psrname=$first
      rajd=`echo $line | awk '{print $2;}'`
      decjd=`echo $line | awk '{print $3;}'`
      uxt=`echo $line | awk '{print int($4);}'` # integer
      uxtime=`echo $line | awk '{print $4;}'`   # float
      freq=`echo $line | awk '{print $5;}'`
      
      sens_file=${psrname}_${rajd}_${decjd}_${uxt}_${freq}MHz_sensitivity.out
      
      echo "python ~/github/station_beam/python/eda_sensitivity.py --bandwidth $bw --frequency $freq --unixtime $uxt --ra=$rajd --dec=$decjd --inttime=$inttime --outsens_file=uxtime${uxt}_eda2_sensitivity --antenna_locations=${HOME}/aavs-calibration/config/eda2/antenna_locations.txt --outfile_mode=a  --trcv_type=trcv_from_skymodel_with_err --nos11 --header=HEADER --use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2 -p None -m analytic ${options} > ${sens_file} 2>&1"
      python ~/github/station_beam/python/eda_sensitivity.py \
         --bandwidth $bw \
         --frequency $freq \
         --unixtime $uxt \
         --ra=$rajd \
         --dec=$decjd \
         --inttime=$inttime \
         --outsens_file=uxtime${uxt}_eda2_sensitivity \
         --antenna_locations=${HOME}/aavs-calibration/config/eda2/antenna_locations.txt \
         --outfile_mode=a \
         --trcv_type=trcv_from_skymodel_with_err \
         --nos11 \
         --header=HEADER \
         --use_beam_fits \
         --station_name=EDA \
         --size=512 \
         --trcv_type=trcv_eda2 \
         -p None \
         -m analytic ${options} \
         > ${sens_file} 2>&1

    sefd_x=`tail -6 ${sens_file} | grep SEFD_XX | awk '{print $25;}'`
    sefd_y=`tail -6 ${sens_file} | grep SEFD_YY | awk '{print $25;}'`
    sefd_i=`echo $sefd_x" "$sefd_y | awk '{print 0.5*sqrt($1*$1+$2*$2);}'`
    
    outline="$line $sefd_x $sefd_y $sefd_i"
    echo $outline >> $outfile
   else
#      echo "INFO : comment line |$line| skipped"
      outline="$line SEFD_X SEFD_Y SEFD_I(=0.5*sqrt(SEFD_X^2 + SEFD_Y^2))"
      echo $outline >> $outfile
   fi
done < $infile