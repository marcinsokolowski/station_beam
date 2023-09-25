#!/bin/bash

# BASED ON : 
#   /data4/eda/eda2/2019_06_08_drift_scan_zenith/
#   

start_time=20190626_160000
start_time_ux=`date2date -local2ux=${start_time}`
if [[ -n "$1" && "$1" != "-" ]]; then
   start_time=$1
   len=`echo $start_time | awk '{print length($1);}'`
   is_unixtime=`echo $start_time | awk '{if(length($1)==10 && int($1)>1000000000){print 1;}else{print 0;}}'`

   if [[ $is_unixtime -gt 0 ]]; then
      start_time_ux=$start_time
      start_time=`date -d "1970-01-01 UTC $start_time_ux seconds" +"%Y%m%d_%T"`
   else   
      if [[ $len -ge 15 ]]; then
         start_time_ux=`date2date -local2ux=${start_time}`
      else
         if [[ $start_time -lt 20500101 ]]; then
            start_time_ux=`date2date -local2ux=${start_time}_160000`
         fi
      fi
  fi
fi

echo "$start_time -> $start_time_ux"

n_days=30
if [[ -n "$2" && "$2" != "-" ]]; then
   n_days=$2
fi

timeres=60
if [[ -n "$3" && "$3" != "-" ]]; then
   timeres=$3
fi

simul_dir=/raid/data/eda/eda2/simulations
if [[ -n "$4" && "$4" != "-" ]]; then
   simul_dir=$4
fi

FEKO_pattern_file=EDA2_48elem_station_beam160MHz_Xpol.txt # ideally in /data4/eda/eda2/simulation/AntennaPatterns/
if [[ -n "$5" && "$5" != "-" ]]; then
   FEKO_pattern_file=$5
fi

antenna_location_file=$BIGHORNS/software/analysis/scripts/shell/eda/eda2/config/antenna_locations_eda2_ph1.txt
if [[ -n "$6" && "$6" != "-" ]]; then
   antenna_location_file=$6
fi

parallel=0
if [[ -n "$7" && "$7" != "-" ]]; then
   parallel=$7
fi

if [[ $parallel -gt 0 ]]; then
   if [[ $n_days -gt 1 ]]; then
      echo "WARNING : n_days = $n_days > 1 -> cannot use parallel mode = 1 -> chaning back to 0"
      parallel=0
   fi
fi

freq_mhz=159
if [[ -n "$8" && "$8" != "-" ]]; then
   freq_mhz=$8
fi
freq_mhz_str=`echo $freq_mhz | awk '{printf("%03d\n",$1);}'`
freq_cc=`echo $freq_mhz | awk '{printf("%d\n",$1/1.28);}'`

rms_phase_error=-1000
rms_error_option=""
if [[ -n "$9" && "$9" != "-" ]]; then
   rms_phase_error=$9
   rms_error_option="--gain_sigma_ph=${rms_phase_error}"
fi


if [[ -n "${10}" && "${10}" != "-" ]]; then
   rms_error_option="${10}"
fi


pointing_az_deg=0.00
if [[ -n "${11}" && "${11}" != "-" ]]; then
   pointing_az_deg=${11}
fi

pointing_za_deg=0.00
if [[ -n "${12}" && "${12}" != "-" ]]; then
   pointing_za_deg=${12}
fi

do_sensitivity=1
if [[ -n "${13}" && "${13}" != "-" ]]; then
   do_sensitivity=${13}
fi

station_name=EDA2
if [[ -n "${14}" && "${14}" != "-" ]]; then
   station_name=${14}
fi
station_name_lower=`echo $station_name | awk '{print tolower($1);}'`

model_options=""
if [[ -n "${15}" && "${15}" != "-" ]]; then
   model_options=${15}
fi


echo "#################################################################################"
echo "PARAMETERS :"
echo "#################################################################################"
echo "station_name = $station_name ( $station_name_lower )"
echo "start_time = $start_time ( start_time_ux = $start_time_ux )"
echo "n_days     = $n_days"
echo "timeres    = $timeres"
echo "simul_dir  = $simul_dir"
echo "( pointing_az_deg , pointing_za_deg ) = ( $pointing_az_deg , $pointing_za_deg ) [deg]"
echo "FEKO_pattern_file     = $FEKO_pattern_file"
echo "antenna_location_file = $antenna_location_file"
echo "parallel              = $parallel"
echo "frequency   = $freq_mhz [MHz] ( string version = $freq_mhz_str) "
echo "rms error   = $rms_phase_error [deg] (<0 - means not used) -> option = $rms_error_option"
echo "sensitivity = ${do_sensitivity}"
echo "model_options = $model_options"
echo "#################################################################################"

cd $simul_dir
pwd

day=0
uxtime=$start_time_ux
while [[ $day -le $n_days ]];
do
   start_time=`date -d "1970-01-01 UTC $uxtime seconds" +"%Y%m%d_%H%M%S"`
   echo
   echo "Generating day - $day , uxtime = $uxtime , start_time = $start_time"

   mkdir -p ${start_time}
   cd ${start_time}
   pwd
   echo "gentotal_new! $start_time $timeres $n_days $freq_mhz"
   gentotal_new! $start_time $timeres $n_days $freq_mhz

   # BASED ON : /data4/eda/eda2/2019_06_08_drift_scan_zenith/

   # HASLAM MAP and simple analytic model for the EDA-2 :   
   if [[ $parallel -gt 0 ]]; then
       echo "nohup simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 \"-q -l -b\" - - - \"--plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $rms_error_option $model_options\" > ${start_time}_analytic.out 2>&1 &"
       nohup simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 "-q -l -b" - - - "--plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $rms_error_option $model_options" > ${start_time}_analytic.out 2>&1 &

       sleep 2
   else
       echo "simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 \"-q -l -b\" - - - \"--plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $rms_error_option $model_options\" > ${start_time}_analytic.out 2>&1"
       simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 "-q -l -b" - - - "--plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $rms_error_option $model_options" > ${start_time}_analytic.out 2>&1 
   fi   

   #  Showspec/GSM + FEKO from Daniel : /data4/eda/eda2/simulation/AntennaPatterns/EDA2_48elem_station_beam160MHz_Xpol.txt
   # based on : bighorns@bighorns:/data4/eda/eda2/2019_06_08_drift_scan_zenith/showspec
   ls -al ${FEKO_pattern_file}
   if [[ -s ${FEKO_pattern_file} || -e ${FEKO_pattern_file} ]]; then
       if [[ $rms_phase_error -le 0 ]]; then
           ln -s /data/aavs/data/20190411/simulation/gsm   
           ln -s ${FEKO_pattern_file}
           if [[ $parallel -gt 0 ]]; then
               echo "nohup showspec_aavs1! ${uxtime} ${FEKO_pattern_file} 1 $timeres 86400 $freq_mhz > showspec_aavs1.out 2>&1 &"
               nohup showspec_aavs1! ${uxtime} ${FEKO_pattern_file} 1 $timeres 86400 $freq_mhz > showspec_aavs1.out 2>&1 &
               sleep 2
           else
               echo "showspec_aavs1! ${uxtime} ${FEKO_pattern_file} 1 $timeres 86400 $freq_mhz > showspec_aavs1.out 2>&1"
               showspec_aavs1! ${uxtime} ${FEKO_pattern_file} 1 $timeres 86400 $freq_mhz > showspec_aavs1.out 2>&1
           fi       
   
           # HASLAM MAP and FEKO model :
           # based on : bighorns@bighorns:/data4/eda/eda2/2019_06_08_drift_scan_zenith/FEKO_HASLAM
           mkdir FEKO_HASLAM
           cd FEKO_HASLAM
           echo "ln -s ../total_power_${freq_mhz_str}_${freq_mhz_str}_MHz.txt"
           ln -s ../total_power_${freq_mhz_str}_${freq_mhz_str}_MHz.txt
           ln -s ${FEKO_pattern_file}
   
           if [[ $parallel -gt 0 ]]; then
               echo "nohup simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 \"-q -l -b\" - - - \"--size=5000 --plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $model_options\" - - - - - ${FEKO_pattern_file} > simulate_total_power_vs_time_all 2>&1 &"
               nohup simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 "-q -l -b" - - - "--size=5000 --plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $model_options" - - - - - ${FEKO_pattern_file} > simulate_total_power_vs_time_all 2>&1 &
           else
               echo "simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 \"-q -l -b\" - - - \"--size=5000 --plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $model_options\" - - - - - ${FEKO_pattern_file} > simulate_total_power_vs_time_all 2>&1"
               simulate_total_power_vs_time_all.sh $timeres $pointing_az_deg $pointing_za_deg 1 "-q -l -b" - - - "--size=5000 --plottype=none --antenna_list_file=$antenna_location_file --antlist_start_column=1 $model_options" - - - - - ${FEKO_pattern_file} > simulate_total_power_vs_time_all 2>&1 
           fi
           cd ../
       else
           echo "WARNING : phase error is not implemented when FEKO file is provided (does not make sense in this case)"
       fi
   else
       echo "WARNING : the specified FEKO file does not exist -> simulation of lightcurves using FEKO files skipped"
   fi
   
   mkdir sensitivity/
   cd sensitivity/
   if [[ $do_sensitivity -gt 0 ]]; then 
      export PATH=~/github/station_beam/:$PATH
   
      station_name2=${station_name}
      if [[ $station_name == "EDA2" ]]; then
         station_name2="EDA"
      fi

      if [[ $parallel -gt 0 ]]; then
         echo "nohup calc_eda_sensitvity_all.sh $uxtime 86400 - \"--station_name=${station_name2} --size=512 --trcv_type=trcv_${station_name_lower}\" $freq_cc - $timeres 400 400 > out 2>&1 &"
         nohup calc_eda_sensitvity_all.sh $uxtime 86400 - "--station_name=${station_name2} --size=512 --trcv_type=trcv_${station_name_lower}" $freq_cc - $timeres 400 400 > out 2>&1 &
      else
         echo "calc_eda_sensitvity_all.sh $uxtime 86400 - \"--station_name=${station_name2} --size=512 --trcv_type=trcv_${station_name_lower}\" $freq_cc - $timeres 400 400 > out 2>&1"
         calc_eda_sensitvity_all.sh $uxtime 86400 - "--station_name=${station_name2} --size=512 --trcv_type=trcv_${station_name_lower}" $freq_cc - $timeres 400 400 > out 2>&1
      fi
   fi
   cd ../
   
   cd ../           
   uxtime=$(($uxtime+86400))
   day=$(($day+1))
      
   if [[ $parallel -gt 0 ]]; then
      day=$(($n_days+1))
      echo "Only one day run is allowed for the parallel mode -> exiting now (by setting day = $day)"
   fi
done

   