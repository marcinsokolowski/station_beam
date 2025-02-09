#!/bin/bash

# obs_script=doitnoweda2
logfile=J0534+2200_40channels_1800sec_flagants_sepdada_ch256.out
if [[ -n "$1" && "$1" != "-" ]]; then
   logfile="$1"
fi

timestep=300
if [[ -n "$2" && "$2" != "-" ]]; then
   timestep="$2"
fi

do_start_simul=1
if [[ -n "$3" && "$3" != "-" ]]; then
   do_start_simul=$3
fi


# if [[ ! -s ${logfile} ]]; then
#   echo "ERROR: log file $logfile not found -> exiting"
#   exit -1;
#fi

freq_ch_start=256 # `head -20 ${logfile} | grep freq_list | awk '{print $3;}'`
n_channels=40 # `head -20 ${logfile} | grep "N channels" | awk '{print $4;}'`
bw_mhz=`echo $n_channels | awk '{print $1*(400./512);}'`
bw_ch_mhz=`echo 1 | awk '{print $1*(400./512);}'`
freq_ch_end=$(($freq_ch_start+$n_channels))
start_uxtime=1732099753 # `head -20 ${logfile} | grep start_uxtime | awk '{print $3;}'`
interval=3600 # `head -20 ${logfile} | grep interval | awk '{print $3;}'`
end_uxtime=$(($start_uxtime+$interval))
object=FRB20240114  # `head -20 ${logfile} | grep object | awk '{print $3;}'`
rajd=321.9125 # `head -20 ${logfile} | grep "ra       =" | awk '{print $3;}'`
decjd=4.32916667 # `head -20 ${logfile} | grep "dec      =" | awk '{print $3;}'`

# sefd_list_file=${object}_startux${start_uxtime}_${interval}sec_startfreq${freq_ch_start}_${n_channels}chan.txt

echo "-----------------------------------"
echo "OBSERVATION PARAMETERS:"
echo "-----------------------------------"
echo "freq_ch_start = $freq_ch_start"
echo "n_channels    = $n_channels"
echo "bw_mhz        = $bw_mhz"
echo "bw_ch_mhz     = $bw_ch_mhz"
echo "freq_ch_end   = $freq_ch_end"
echo "start_uxtime  = $start_uxtime"
echo "end_uxtime    = $end_uxtime"
echo "interval      = $interval"
echo "object        = $object"
echo "(RA,DEC)      = ($rajd,$decjd) [deg]"
echo "output file   = $sefd_list_file"
echo "timestep      = $timestep [sec]"
echo "do_start_simul = $do_start_simul"
echo "-----------------------------------"

#if [[ -s ${sefd_list_file} ]]; then
#   echo "mv ${sefd_list_file} ${sefd_list_file}.old"
#   mv ${sefd_list_file} ${sefd_list_file}.old
#fi

mv doit! doit.old!

# echo "# OBJECT RA[deg] DEC[deg] UXTIME FREQ[MHz]" > ${sefd_list_file}

uxtime=${start_uxtime}
while [[ $uxtime -le $end_uxtime ]]; 
do
   freq_ch=$freq_ch_start   
   sefd_list_file=${object}_startux${uxtime}_${timestep}sec_startfreq${freq_ch_start}_${n_channels}chan.txt
   echo "# OBJECT RA[deg] DEC[deg] UXTIME FREQ[MHz]" > ${sefd_list_file}
         
   while [[ $freq_ch -le $freq_ch_end ]];
   do
      freq_mhz=`echo $freq_ch | awk '{print $1*(400.00/512.0);}'`
      echo "$object $rajd $decjd $uxtime $freq_mhz" >> ${sefd_list_file}
   
      freq_ch=$(($freq_ch+1))
   done
   
   echo "INFO : generated file $sefd_list_file"
   
   if [[ $do_start_simul -gt 0 ]]; then
      echo "Starting simulations in 10 seconds ..."
      echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}"
      ~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}
      
      echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}" >> doit!
   else
      echo "WARNING : running simulations is not requested"
      echo "In order to do it manually execute the following command (saved in doit! script):"
      echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}"
      echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}" >> doit!
   fi
   
   uxtime=$(($uxtime+$timestep))
done

#if [[ $do_start_simul -gt 0 ]]; then
#  echo "Starting simulations in 10 seconds ..."
#  sleep 10
#  
#   echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}"
#   ~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}
#else
#   echo "WARNING : running simulations is not requested"
#   echo "In order to do it manually execute the following command:"
#   echo "~/github/station_beam/scripts/psrflux/calculate_sefd_for_pointings.sh ${sefd_list_file} - ${timestep} ${bw_ch_mhz}"
#fi








