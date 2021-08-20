#!/bin/bash

# INPUT :
# obs_freq_list.txt = obsID FREQ_cc 

# export PATH="/opt/caastro/ext/anaconda/bin:$PATH"

start_ux=1578441600
if [[ -s last_processed_ux.txt ]]; then
   start_ux=`cat last_processed_ux.txt`
fi
if [[ -n "$1" && "$1" != "-" ]]; then
   start_ux=$1
fi

interval=86400 # 24 hours 
if [[ -n "$2" && "$2" != "-" ]]; then
   interval=$2
fi

out_basename=eda_sensitivity
if [[ -n "$3" && "$3" != "-" ]]; then
   out_basename=$3
fi

# options="--trcv_budi"
options=""
if [[ -n "$4" && "$4" != "-" ]]; then
   options="$4"
fi

freq_cc_list="39 46 54 62 69 70 78 85 93 101 109 117 121 125 132 140 145 148 156 164 169 171 179 187 195 203 210 218 226 234 242 250 257 265 273"
if [[ -n "$5" && "$5" != "-" ]]; then
   freq_cc_list=$5
fi
freq_cc_count=`echo $freq_cc_list | awk '{print $NF;}'`

use_beamf_errors=0
beamf_err_options=""
# beamf_err_options="--gain_sigma_db=0.7 --gain_sigma_ph=8"
if [[ -n "$6" && "$6" != "-" ]]; then
   use_beamf_errors=1
   beamf_err_options="$6"
fi

time_step=60 # seconds
if [[ -n "$7" && "$7" != "-" ]]; then
   time_step=$7
fi

az_step=5 # 1 gave no change !
if [[ -n "$8" && "$8" != "-" ]]; then
   az_step=$8
fi

za_step=5
if [[ -n "$9" && "$9" != "-" ]]; then
   za_step=$9
fi

max_az=360
if [[ -n "${10}" && "${10}" != "-" ]]; then
   max_az=${10}
fi

za_start=0
if [[ -n "${11}" && "${11}" != "-" ]]; then
   za_start=${11}
fi

za_end=90
if [[ -n "${12}" && "${12}" != "-" ]]; then
   za_end=${12}
fi

az_start=0
if [[ -n "${13}" && "${13}" != "-" ]]; then
   az_start=${13}
fi

echo "########################################"
echo "PARAMETERS:"
echo "########################################"
echo "za range : $za_start - $za_end [deg]"
echo "az range : $az_start - $max_az [deg]"
echo "########################################"


eda_sensitivity_path=`which eda_sensitivity.py`

# rm -f ${out_basename}_XX.txt ${out_basename}_YY.txt

end_ux=$(($start_ux+$interval))

ux=$start_ux
while [[ $ux -le $end_ux ]];
do
   gps=`ux2gps! $ux`
   lst=`ux2sid $ux | awk '{print $8;}'` 
   dtm=`ux2ut! $ux`
   echo "------------------ ux = $ux -> gps = $gps -> lst = $lst [h] -> $dtm UTC ------------------"
   obs_list=obs_list_${gps}.txt

   za=${za_start}
   while [[ $za -le ${za_end} ]]; 
   do
      az=${az_start}
      while [[ $az -lt ${max_az} ]]; 
      do
          echo "   (az,za) = ($az,$za) [deg]"          
          index=0
          
          if [[ -s ${gps}_az${az}_za${za}_${out_basename}_XX.txt && -s ${gps}_az${az}_za${za}_${out_basename}_YY.txt ]]; then
             echo "Files ${gps}_az${az}_za${za}_${out_basename}_XX.txt and ${gps}_az${az}_za${za}_${out_basename}_YY.txt already exist -> skipped"
          else                                        
             echo "Processing all frequencies for (az,za) = ($az,$za) [deg] , started at :"
             date
 
             ux_start1=`date +%s`         
             for freq_cc in `echo $freq_cc_list`
             do
                 echo "        freq_cc = $freq_cc"
                 header_options=""
                 if [[ $index -gt 0 ]]; then
                    header_options="--no-header"
                 fi              

                 ux_start=`date +%s`
                 echo "python $eda_sensitivity_path -c ${freq_cc} -p None -g ${gps}  -m analytic --az=${az} --za=${za} --outsens_file=${gps}_az${az}_za${za}_${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 --header=HEADER ${header_options} ${options}"
                 python $eda_sensitivity_path -c ${freq_cc} -p None -g ${gps}  -m analytic --az=${az} --za=${za} --outsens_file=${gps}_az${az}_za${za}_${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 --header=HEADER ${header_options} ${options}
                 ux_end=`date +%s`
                 ux_diff=`echo "$ux_end $ux_start" | awk '{print $1-$2;}'`
                 echo "BENCHMARKING : single frequency step took $ux_diff seconds"
          
                 index=$(($index+1))
             done
             ux_end1=`date +%s`
             ux_diff1=`echo "$ux_end1 $ux_start1" | awk '{print $1-$2;}'`
             echo "BENCHMARKING : simulation of $freq_cc_count frequency channels took $ux_diff1 seconds"

          fi
          az=$(($az+$az_step))
          
          if [[ $za == 0 ]]; then
             echo "INFO : za = $za -> only az=0 calculated -> setting az := $max_az"
             az=${max_az}
          fi
      done   
      za=$(($za+$za_step))
   done
   
   echo $ux > last_processed_ux.txt 
   ux=$(($ux+$time_step))
done

