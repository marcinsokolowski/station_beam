#!/bin/bash

# INPUT :
# obs_freq_list.txt = obsID FREQ_cc 

az_deg=0
if [[ -n "$1" && "$1" != "-" ]]; then
   az_deg=$1
fi

za_deg=0
if [[ -n "$2" && "$2" != "-" ]]; then
   za_deg=$2
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

obs_freq_list_file=obs_freq_list.txt
if [[ -n "$5" && "$5" != "-" ]]; then
   obs_freq_list_file=$5
fi

use_beamf_errors=0
beamf_err_options=""
# beamf_err_options="--gain_sigma_db=0.7 --gain_sigma_ph=8"
if [[ -n "$6" && "$6" != "-" ]]; then
   use_beamf_errors=1
   beamf_err_options="$6"
fi

rm -f ${out_basename}_XX.txt ${out_basename}_YY.txt

eda_sensitivity_path=`which eda_sensitivity.py`

index=0
while read line
do
   obsID=`echo $line | awk '{print $1;}'`
   freq_cc=`echo $line | awk '{print $2;}'`
   
   header_options=""
   if [[ $index -gt 0 ]]; then
      header_options="--no-header"
   fi
   
   echo "python $eda_sensitivity_path -c ${freq_cc} -p None -g ${obsID}  -m analytic --az=${az_deg} --za=${za_deg} --outsens_file=${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 ${header_options} ${options}"
   python $eda_sensitivity_path -c ${freq_cc} -p None -g ${obsID}  -m analytic --az=${az_deg} --za=${za_deg} --outsens_file=${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 ${header_options} ${options}
   
   echo "------------------------------------------------------------ $obsID / $freq_cc ------------------------------------------------------------"
   
   index=$(($index+1))
done < ${obs_freq_list_file}
