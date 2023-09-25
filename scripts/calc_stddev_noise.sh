#!/bin/bash

fp=~/github/station_beam/scripts

station=0 #  0 == AAVS2 and 1 == EDA2
if [[ -n "$1" && "$1" != "-" ]]; then
   station=$1
fi

ar_file=20230107T174258_ch256.ar
if [[ -n "$2" && "$2" != "-" ]]; then
   ar_file=$2
fi

pdv -FTtp ${2} > ascii_archive.txt
# info=$(python ${fp}/get_sigma_clip_manual.py $3 $4)
# std_vcs=$(echo $info | awk '{print $1}')
# off_pts=$(echo $info | awk '{print $2}')
# prof_sum=$(echo $info | awk '{print $3}')
metadata=$(vap -c "name nbin mjd freq period length bw" -n ${2} 2>/dev/null)
src=$(echo $metadata | awk '{print $2}')
bins=$(echo $metadata | awk '{print $3}')
mjd=$(echo $metadata | awk '{print $4}')
freq=$(echo $metadata | awk '{print $5}')
period=$(echo $metadata | awk '{print $6}')
inttime=$(echo $metadata | awk '{print $7}')
bw_MHz=$(echo $metadata | awk '{print $8}')
bw=$(echo "$bw_MHz*1000000" | bc)
uxt=$(python ${fp}/mjd_uxt.py $mjd)
psrinfo=$(psrcat -c "RAJD DECJD" -o short -nohead -nonumber $src)
rajd=$(echo $psrinfo | awk '{print $1}')
decjd=$(echo $psrinfo | awk '{print $2}')


if [[ $station == 0 ]]; then 
   # AAVS2 :
   echo "python ~/github/station_beam//eda_sensitivity.py --bandwidth $bw --frequency $freq --unixtime $uxt --ra=$rajd --dec=$decjd --inttime=$inttime --outsens_file=uxtime${uxt}_aavs2_sensitivity --antenna_locations=/home/msok/aavs-calibration/config/aavs2/antenna_locations.txt --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err --nos11 --header=HEADER --use_beam_fits --station_name=SKALA4 --size=512 --trcv_type=trcv_aavs2_vs_za_deg -p None -m analytic"
							
   python ~/github/station_beam/python/eda_sensitivity.py \
					--bandwidth $bw \
					--frequency $freq \
					--unixtime $uxt \
					--ra=$rajd \
					--dec=$decjd \
					--inttime=$inttime \
					--outsens_file=uxtime${uxt}_aavs2_sensitivity \
					--antenna_locations=/home/msok/aavs-calibration/config/aavs2/antenna_locations.txt \
					--outfile_mode=a \
					--trcv_type=trcv_from_skymodel_with_err \
					--nos11 \
					--header=HEADER \
					--use_beam_fits \
					--station_name=SKALA4 \
					--size=512 \
					--trcv_type=trcv_aavs2_vs_za_deg \
					-p None \
					-m analytic \
					> ${src}_sensitivity.out 2>&1
else					
   echo "python ~/github/station_beam//eda_sensitivity.py --bandwidth $bw --frequency $freq --unixtime $uxt --ra=$rajd --dec=$decjd --inttime=$inttime --outsens_file=uxtime${uxt}_eda2_sensitivity --antenna_locations=/home/msok/aavs-calibration/config/eda2/antenna_locations.txt --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err --nos11 --header=HEADER --use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2 -p None -m analytic"

   python ~/github/station_beam/python/eda_sensitivity.py \
					--bandwidth $bw \
					--frequency $freq \
					--unixtime $uxt \
					--ra=$rajd \
					--dec=$decjd \
					--inttime=$inttime \
					--outsens_file=uxtime${uxt}_eda2_sensitivity \
					--antenna_locations=/home/msok/aavs-calibration/config/eda2/antenna_locations.txt \
					--outfile_mode=a \
					--trcv_type=trcv_from_skymodel_with_err \
					--nos11 \
					--header=HEADER \
					--use_beam_fits \
					--station_name=EDA \
					--size=512 \
					--trcv_type=trcv_eda2 \
					-p None \
					-m analytic \
					> ${src}_sensitivity.out 2>&1
fi					

lastline=$(grep 'Stokes I images' ${src}_sensitivity.out)
std_sim=$(echo $lastline | awk '{print $8}')

echo "Output of simulation in file : ${src}_sensitivity.out"
echo $lastline