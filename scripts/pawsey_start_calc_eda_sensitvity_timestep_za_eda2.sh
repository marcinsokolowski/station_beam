#!/bin/bash -l

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
module load astropy/v2.0.12

echo "Modules avail:"
module avail

echo
echo "Modules list:"
module list

start_ux=$1
za_deg=$2

param_file="parameters.txt"
station=EDA2
path="/astro/mwaops/msok/eda2/simulations/${station}"

mkdir -p ${path}
cd ${path}

# echo "start_ux = $start_ux" > ${param_file}
# echo "out_basename = eda_sensitivity" >> ${param_file}
# echo "options = --use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" >> ${param_file}
# echo "freq_cc_list = 39 46 54 62 69 70 78 85 93 101 109 117 121 125 132 140 145 148 156 164 169 171 179 187 195 203 210 218 226 234 242 250 257 265 273"  >> ${param_file}
# echo "za_param = $za_deg" >> ${param_file}
# echo "az_step = 5" >> ${param_file} 
# echo "za_step = 5" >> ${param_file} 
# echo "max_az  = 360" >> ${param_file}



echo "Started at :" > ${start_ux}.log
date >> ${start_ux}.log

echo "Generated parameter file ${param_file} :"
echo "----------------------------------------"
cat ${param_file}
echo "----------------------------------------"

pwd
echo "$run calc_eda_sensitvity_timestep_za.sh $start_ux eda_sensitivity \"--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2\" > ${start_ux}.out"
$run calc_eda_sensitvity_timestep_za.sh $start_ux eda_sensitivity "--use_beam_fits --station_name=EDA --size=512 --trcv_type=trcv_eda2" > ${start_ux}.out 

