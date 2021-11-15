#!/bin/bash

object_radec()
{
   object=$1
   
   ret=""
   if [[ $object == "3c444" || $object == "3C444" ]]; then
      ra=333.607249950
      dec=-17.02661111
   fi 
 
   if [[ $object == "eor0" || $object == "EOR0" || $object == "EoR0" ]]; then  
      ra=0
      dec=-27
   fi
   
   if [[ $object == "eor1" || $object == "EOR1" || $object == "EoR1" ]]; then  
      ra=60
      dec=-27
   fi

   if [[ $object == "hyda" || $object == "HYDA" || $object == "HYDRAA" || $object == "HYDRA-A" || $object == "hydra-a" ]]; then  
      ra=139.52374995
      dec=-12.09555556
   fi

   if [[ $object == "vira" || $object == "VIRA" || $object == "VIRGOA" || $object == "VIRGO-A" || $object == "virgo-a" ]]; then  
      ra=187.70416667
      dec=12.39111111
   fi
   
   echo "DEBUG : object_radec -> object = $object, (ra,dec) = ($ra,$dec) [deg]"
}   


object=3c444
if [[ -n "$1" && "$1" != "-" ]]; then
   object=$1
fi

freq_mhz=160
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_mhz=$2
fi

gps_start=1320812567
if [[ -n "$3" && "$3" != "-" ]]; then
   gps_start=$3
fi

# setup RA,DEC by the source name as default :
ra=333.607249950 # degree
dec=-17.02661111 # degree
object_radec $object

# Can overwrite RA,DEC with parameters:
if [[ -n "$4" && "$4" != "-" ]]; then
   ra=$4
fi

# Can overwrite RA,DEC with parameters:
if [[ -n "$5" && "$5" != "-" ]]; then
   dec=$5
fi

# Hardcoded (for now)
duration=86400
gps_end=$(($gps_start+$duration))
step=300
station=aavs2
staton_name=SKALA4
if [[ $station == "eda2" ]]; then
   staton_name=EDA2
fi

echo "############################################################"
echo "PARAMETERS :"
echo "############################################################"
echo "object = $object"
echo "(RA,DEC) = ($ra,$dec) [degree]"
echo "Frequency = $freq [MHz]"
echo "GPS START = $gps_start"
echo "Duration  = $duration [seconds]"
echo "Step      = $step [seconds]"
echo "############################################################"

if [[ ! -s antenna_locations_${station}.txt ]]; then
   echo "cp ~/aavs-calibration/config/${station}/antenna_locations.txt antenna_locations_${station}.txt"
   cp ~/aavs-calibration/config/${station}/antenna_locations.txt antenna_locations_${station}.txt
else
   echo "Antenna location file antenna_locations_${station}.txt already exists"   
fi

gps=$gps_start
while [[ $gps -le $gps_end ]];
do
   echo "python ~/github/station_beam//eda_sensitivity.py --freq=${freq_mhz} -p None -g ${gps}  -m analytic --ra=${ra} --dec=${dec} --outsens_file=${object}_${station}_sensitivity --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err  --nos11 --header=HEADER  --use_beam_fits --station_name=${staton_name} --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_${station}.txt --projection=aee"
   python ~/github/station_beam//eda_sensitivity.py --freq=${freq_mhz} -p None -g ${gps}  -m analytic --ra=${ra} --dec=${dec} --outsens_file=${object}_${station}_sensitivity --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err  --nos11 --header=HEADER  --use_beam_fits --station_name=${staton_name} --size=512 --trcv_type=trcv_aavs2_vs_za_deg --antenna_locations=antenna_locations_${station}.txt --projection=aee

   gps=$(($gps+$step))
done   

# elev vs. azimuth :
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){azim=$10;elev=(90.00-$15);}}else{if(azim>180){azim=azim-360;}if(azim>=-120 && azim<=120){print azim" "elev;}}}' ${object}_${station}_sensitivity_XX.txt > ${object}_elev_vs_azimuth.txt

# get A/T vs. azimuth 
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$2;}}}' ${object}_${station}_sensitivity_XX.txt | sort -n > ${object}_${station}_sensitivity_vs_azim_XX.txt
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$2;}}}' ${object}_${station}_sensitivity_YY.txt | sort -n > ${object}_${station}_sensitivity_vs_azim_YY.txt
ls ${object}_${station}_sensitivity_vs_azim_XX.txt ${object}_${station}_sensitivity_vs_azim_YY.txt > ${object}_${station}_sensitivity_vs_azim.list

# get A/T vs. elevation :
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=90.00-$15;}}else{if(elev>0){print elev" "$2;}}}' ${object}_${station}_sensitivity_XX.txt > ${object}_${station}_sensitivity_vs_elev_XX.txt
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=90.00-$15;}}else{if(elev>0){print elev" "$2;}}}' ${object}_${station}_sensitivity_YY.txt > ${object}_${station}_sensitivity_vs_elev_YY.txt
ls ${object}_${station}_sensitivity_vs_elev_XX.txt ${object}_${station}_sensitivity_vs_elev_YY.txt > ${object}_${station}_sensitivity_vs_elev.list

# get A vs. azimuth 
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$4;}}}' ${object}_${station}_sensitivity_XX.txt | sort -n > ${object}_eda_aeff_vs_azim_XX.txt
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$4;}}}' ${object}_${station}_sensitivity_YY.txt | sort -n > ${object}_eda_aeff_vs_azim_YY.txt
ls ${object}_eda_aeff_vs_azim_XX.txt ${object}_eda_aeff_vs_azim_YY.txt > ${object}_eda_aeff_vs_azim.list

# get T_sys vs. azimuth 
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$3;}}}' ${object}_${station}_sensitivity_XX.txt | sort -n > ${object}_eda_tsys_vs_azim_XX.txt
awk -v elev=-1 '{if($1=="#"){if($2=="HEADER"){elev=$10;}}else{if(elev>180){elev=elev-360;}if(elev>=-120 && elev<=120){print elev" "$3;}}}' ${object}_${station}_sensitivity_YY.txt | sort -n > ${object}_eda_tsys_vs_azim_YY.txt
ls ${object}_eda_tsys_vs_azim_XX.txt ${object}_eda_tsys_vs_azim_YY.txt > ${object}_eda_tsys_vs_azim.list

# plots may not work without root package installed :
mkdir -p images/
root -b -q -l "plotNfiles_AoT_vs_azim_PAPER_LFAAREQ_TEST.C(\"${object}_${station}_sensitivity_vs_azim.list\",${freq_mhz},\"${object}_elev_vs_azimuth.txt\")"
root -b -q -l "plotNfiles_AoT_vs_elev_PAPER_LFAAREQ_TEST.C(\"${object}_${station}_sensitivity_vs_elev.list\",${freq_mhz})"

root -b -q -l "plotNfiles_Aeff_vs_azim_PAPER_LFAAREQ_TEST.C(\"${object}_eda_aeff_vs_azim.list\",${freq_mhz})"
root -b -q -l "plotNfiles_Tsys_vs_azim_PAPER_LFAAREQ_TEST.C(\"${object}_eda_tsys_vs_azim.list\",${freq_mhz})"

# checks :
awk -v cnt=0 -v mean=0.00 -v elev=-1 -v azim=-1 '{if($1=="#"){if($2=="HEADER"){azim=$10;elev=(90.00-$15);}}else{if(azim>180){azim=azim-360;}if(elev>=45){mean+=$2;cnt+=1;print azim" "elev" "$2" "mean/cnt;}}}' ${object}_${station}_sensitivity_XX.txt


