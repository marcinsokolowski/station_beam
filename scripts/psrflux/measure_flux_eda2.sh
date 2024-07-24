#!/bin/bash

# Based on script developed by Christopher Lee (paper : https://ui.adsabs.harvard.edu/abs/2022PASA...39...42L/abstract)

fp=~/github/station_beam/scripts/psrflux/
# export PATH=~/github/station_beam/scripts/psrflux/:$PATH

ar_file=20240501T101103_ch290.ar
if [[ -n "$1" && "$1" != "-" ]]; then
   ar_file="$1"
fi

pulse_start_bin=0
if [[ -n "$2" && "$2" != "-" ]]; then
   pulse_start_bin=$2
fi

pulse_end_bin=0
if [[ -n "$3" && "$3" != "-" ]]; then
   pulse_end_bin=$3
fi

echo "#########################################"
echo "PARAMETERS:"
echo "#########################################"
echo "ar_file = $ar_file"
echo "pulse_start_bin = $pulse_start_bin"
echo "pulse_end_bin = $pulse_end_bin"
echo "#########################################"


pdv -FTtp ${ar_file} > ascii_archive.txt

echo "python $fp/get_sigma_clip_manual.py $pulse_start_bin $pulse_end_bin"
info=$(python $fp/get_sigma_clip_manual.py $pulse_start_bin $pulse_end_bin)
std_vcs=$(echo $info | awk '{print $1}')
off_pts=$(echo $info | awk '{print $2}')
prof_sum=$(echo $info | awk '{print $3}')

metadata=$(vap -c "name nbin mjd freq period length bw" -n ${ar_file} 2>/dev/null)
src=$(echo $metadata | awk '{print $2}')
bins=$(echo $metadata | awk '{print $3}')
mjd=$(echo $metadata | awk '{print $4}')
freq=$(echo $metadata | awk '{print $5}')
period=$(echo $metadata | awk '{print $6}')
inttime=$(echo $metadata | awk '{print $7}')
bw_MHz=$(echo $metadata | awk '{print $8}')
bw=$(echo "$bw_MHz*1000000" | bc)

pulse_profile_file=${src}_${ar_file%%.ar}_pulse_profile.txt
sens_file=${src}_${ar_file%%.ar}_sensitivity.out
awk '{if(NR>1){print $3" "$4;}}' ascii_archive.txt > ${pulse_profile_file}

echo ""python $fp/mjd_uxt.py $mjd
uxt=$(python $fp/mjd_uxt.py $mjd)

psrinfo=$(psrcat -c "RAJD DECJD" -o short -nohead -nonumber $src)
rajd=$(echo $psrinfo | awk '{print $1}')
decjd=$(echo $psrinfo | awk '{print $2}')

echo "python ~/github/station_beam/python/eda_sensitivity.py --bandwidth $bw --frequency $freq --unixtime $uxt --ra=$rajd --dec=$decjd --inttime=$inttime --outsens_file=uxtime${uxt}_eda2_sensitivity --antenna_locations=${HOME}/aavs-calibration/config/eda2/antenna_locations.txt --outfile_mode=a  --trcv_type=trcv_from_skymodel_with_err --nos11 --header=HEADER --use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2 -p None -m analytic > ${sens_file} 2>&1"
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
         -m analytic \
         > ${sens_file} 2>&1

lastline=$(grep 'Stokes I images' ${sens_file})
std_sim=$(echo $lastline | awk '{print $8}')

#Debugging
echo "Source: $src"
echo "Centre frequency: $freq"
echo "Integration time: $inttime"
echo "Middle unix time: $uxt"
echo "MJD: $mjd"
echo "Right ascention: $rajd"
echo "Declination: $decjd"
echo "Period: $period"

echo "PARAMETERS of flux_calc.py:"
echo "std_vcs = $std_vcs"
echo "std_sim = $std_sim"
echo "bins    = $bins"
echo "off_pts = $off_pts"
echo "period  = $period"
echo "prof_sum = $prof_sum"

echo "python ${fp}/flux_calc.py $std_vcs $std_sim $bins $off_pts $period $prof_sum"
python ${fp}/flux_calc.py $std_vcs $std_sim $bins $off_pts $period $prof_sum


echo "MSOK's calculation in root"
echo "plot_psr_profile.C(\"${pulse_profile_file}\",true,\"pulse_gauss\",0,1.0,$std_sim)"
echo "WARNING : modify noise range manually !!!"
cp ~/github/station_beam/scripts/psrflux/plot_psr_profile.C .
root -l "plot_psr_profile.C(\"${pulse_profile_file}\",true,\"pulse_gauss\",0,1.0,$std_sim)"
