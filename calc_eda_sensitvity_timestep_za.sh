#!/bin/bash

# script for PAWSEY style , reads parameters from file parameters.txt in format :
# start_ux = 1578441600

# INPUT :
# obs_freq_list.txt = obsID FREQ_cc 

# export PATH="/opt/caastro/ext/anaconda/bin:$PATH"

# last_ux=1578441000
# if [[ -s last_processed_ux.txt ]]; then
#   last_ux=`cat last_processed_ux.txt`
#fi

param_file="parameters.txt"

start_ux=1578441600
out_basename=eda_sensitivity
options=""
freq_cc_list="39 46 54 62 69 70 78 85 93 101 109 117 121 125 132 140 145 148 156 164 169 171 179 187 195 203 210 218 226 234 242 250 257 265 273"
use_beamf_errors=0
beamf_err_options=""
za=0
az_step=5 # 1 gave no change !
max_az=360

pwd
date
#if [[ -s ${param_file} ]]; then
#   start_ux=`cat parameters.txt | grep start_ux | awk '{print $3;}'`
#   out_basename=`cat parameters.txt | grep out_basename | awk '{print $3;}'`
#   options=`cat parameters.txt | grep options | awk 'BEGIN{out="";}{i=index($0,"=");print substr($0,i+1);}'`
#   freq_cc_list=`cat parameters.txt | grep freq_cc_list | awk 'BEGIN{out="";}{i=index($0,"=");print substr($0,i+1);}'`
##   use_beamf_errors=`cat parameters.txt | grep use_beamf_errors | awk 'BEGIN{out="";}{i=index($0,"=");print substr($0,i+1);}'`
##   beamf_err_options=`cat parameters.txt | grep beamf_err_options | awk 'BEGIN{out="";}{i=index($0,"=");print substr($0,i+1);}'`
#   za_param=`cat parameters.txt | grep za_param | awk '{print $3;}'`
#   az_step=`cat parameters.txt | grep az_step | awk '{print $3;}'`
##   za_step=`cat parameters.txt | grep za_step | awk '{print $3;}'`
#   max_az=`cat parameters.txt | grep max_za | awk '{print $3;}'`
#else
#   echo "WARNING : parameter file $param_file not found !"
#   exit;
#fi

if [[ -n "$1" && "$1" != "-" ]]; then
   start_ux=$1
fi

out_basename=eda_sensitivity
if [[ -n "$2" && "$2" != "-" ]]; then
   out_basename=$2
fi

# options="--trcv_budi"
options=""
if [[ -n "$3" && "$3" != "-" ]]; then
   options="$3"
fi

freq_cc_list="39 46 54 62 69 70 78 85 93 101 109 117 121 125 132 140 145 148 156 164 169 171 179 187 195 203 210 218 226 234 242 250 257 265 273"
if [[ -n "$4" && "$4" != "-" ]]; then
   freq_cc_list=$4
fi

# beamf_err_options="--gain_sigma_db=0.7 --gain_sigma_ph=8"
if [[ -n "$5" && "$5" != "-" ]]; then
   use_beamf_errors=1
   beamf_err_options="$5"
fi

if [[ -n "$6" && "$6" != "-" ]]; then
   az_step=$6
fi

if [[ -n "${7}" && "${7}" != "-" ]]; then
   max_az=${7}
fi

python_path=`which python`
# python_path=/pawsey/cle60up07/apps/gcc/4.8.5/python/2.7.14/bin/python

echo "########################################"
echo "PARAMETERS:"
echo "########################################"
echo "start_ux          = $start_ux"
echo "za                = $za_param"
echo "out_basename      = $out_basename"
echo "options           = $options"
echo "freq_cc_list      = $freq_cc_list"
echo "use_beamf_errors  = $use_beamf_errors"
echo "beamf_err_options = $beamf_err_options"
echo "az_step           = $az_step"
echo "max_az            = $max_az"
echo "python_path       = $python_path"
echo "########################################"


eda_sensitivity_path=`which eda_sensitivity.py`

# rm -f ${out_basename}_XX.txt ${out_basename}_YY.txt

gps=`ux2gps! $ux`
lst=`ux2sid $ux | awk '{print $8;}'` 
dtm=`ux2ut! $ux`
echo "------------------ ux = $ux -> gps = $gps -> lst = $lst [h] -> $dtm UTC ------------------"
obs_list=obs_list_${gps}.txt

za=${za_param}

az=0
while [[ $az -lt ${max_az} ]]; 
do
    echo "   (az,za) = ($az,$za) [deg]"          
    index=0
          
    if [[ -s ${gps}_az${az}_za${za}_${out_basename}_XX.txt && -s ${gps}_az${az}_za${za}_${out_basename}_YY.txt ]]; then
       echo "Files ${gps}_az${az}_za${za}_${out_basename}_XX.txt and ${gps}_az${az}_za${za}_${out_basename}_YY.txt already exist -> skipped"
    else                                        
       echo "Processing all frequencies for (az,za) = ($az,$za) [deg] , started at :"
       date
       
       for freq_cc in `echo $freq_cc_list`
       do
           echo "        freq_cc = $freq_cc"
           header_options=""
           if [[ $index -gt 0 ]]; then
              header_options="--no-header"
           fi              
           echo "$python_path $eda_sensitivity_path -c ${freq_cc} -p None -g ${gps}  -m analytic --az=${az} --za=${za} --outsens_file=${gps}_az${az}_za${za}_${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 --header=HEADER ${header_options} ${options}"
           $python_path $eda_sensitivity_path -c ${freq_cc} -p None -g ${gps}  -m analytic --az=${az} --za=${za} --outsens_file=${gps}_az${az}_za${za}_${out_basename} --outfile_mode=a --trcv_type=trcv_from_skymodel_with_err ${beamf_err_options} --nos11 --header=HEADER ${header_options} ${options}
       
           index=$(($index+1))
       done
    fi
    az=$(($az+$az_step))
          
    if [[ $za == 0 ]]; then
       echo "INFO : za = $za -> only az=0 calculated -> setting az := $max_az"
       az=${max_az}
    fi
done   
   
# echo $ux > last_processed_ux.txt 

