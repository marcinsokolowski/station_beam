#!/bin/bash

# Cotter the data
#SBATCH --account=mwaops
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=64gb
#SBATCH --output=/astro/mwaops/msok/log/sensitivity/calc_sensitivity.o%j
#SBATCH --error=/astro/mwaops/msok/log/sensitivity/calc_sensitivity.e%j
#SBATCH --export=NONE

cores=1
run="srun --nodes=1 -c $cores --export=all"

export PATH=/home/msok/github/station_beam:$PATH

# script for PAWSEY style , reads parameters from file parameters.txt in format :
# start_ux = 1578441600

# INPUT :
# obs_freq_list.txt = obsID FREQ_cc 

# export PATH="/opt/caastro/ext/anaconda/bin:$PATH"

# last_ux=1578441000
# if [[ -s last_processed_ux.txt ]]; then
#   last_ux=`cat last_processed_ux.txt`
#fi

# LOAD MODULES :
# module use /group/mwa/software/modulefiles
# module load MWA_Tools/mwa-sci
# module load astropy
# load modules :
module use /group/mwa/software/modulefiles
# module load MWA_Tools/mwa-sci
# module load pawseytools
# module load python/2.7.14
# module load astropy
module load numpy/1.13.3
module load funcsigs/1.0.2
module load setuptools/38.2.1
module load pytest/3.3.0
module load py/1.5.2
module load astropy/v2.0.6
module load ephem
module load scipy
module load matplotlib

# echo "Modules avail:"
# module avail
# echo
# echo "Modules list:"
# module list


param_file="parameters.txt"

start_ux=1578441600
if [[ -n "$1" && "$1" != "-" ]]; then
   start_ux=$1
fi

za=0
if [[ -n "$2" && "$2" != "-" ]]; then
   za=$2
fi

out_basename=eda_sensitivity
if [[ -n "$3" && "$3" != "-" ]]; then
   out_basename=$3
fi

# options="--trcv_budi"
options="--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2"
if [[ -n "$4" && "$4" != "-" ]]; then
   options="$4"
fi

freq_cc_list="39 46 54 62 69 70 78 85 93 101 109 117 121 125 132 140 145 148 156 164 169 171 179 187 195 203 210 218 226 234 242 250 257 265 273"
if [[ -n "$5" && "$5" != "-" ]]; then
   freq_cc_list=$5
fi

use_beamf_errors=0
beamf_err_options=""
# beamf_err_options="--gain_sigma_db=0.7 --gain_sigma_ph=8"
if [[ -n "$6" && "$6" != "-" ]]; then
   use_beamf_errors=1
   beamf_err_options="$6"
fi

az_step=5
if [[ -n "$7" && "$7" != "-" ]]; then
   az_step=$7
fi

max_az=360
if [[ -n "${8}" && "${8}" != "-" ]]; then
   max_az=${8}
fi

python_path=`which python`
# python_path=/pawsey/cle60up07/apps/gcc/4.8.5/python/2.7.14/bin/python

echo "########################################"
echo "PARAMETERS:"
echo "########################################"
echo "start_ux          = $start_ux"
echo "za                = $za"
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

ux=$start_ux
gps=`/home/msok/github/station_beam/ux2gps! $ux`
lst=`/home/msok/mwa_soft/bin/ux2sid $ux | awk '{print $8;}'` 
dtm=`/home/msok/github/station_beam/ux2ut! $ux`
echo "------------------ ux = $ux -> gps = $gps -> lst = $lst [h] -> $dtm UTC ------------------"
obs_list=obs_list_${gps}.txt

za=${za}
echo "za = $za"

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

