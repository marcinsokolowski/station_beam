#!/bin/bash

skymodel=MODEL
if [[ -n "$8" && "$8" != "-" ]]; then
   skymodel=$8
fi

total_power_file=total_power_080_090_MHz.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   total_power_file=$1
fi

outdir=${total_power_file%%.txt}_${skymodel}/
if [[ -n "$2" && "$2" != "-" ]]; then
   outdir=$2
fi

step=300
if [[ -n "$3" && "$3" != "-" ]]; then
   step=$3
fi   

pointing_az_deg=0
if [[ -n "$4" && "$4" != "-" ]]; then
   pointing_az_deg=$4
fi
   
pointing_za_deg=0
if [[ -n "$5" && "$5" != "-" ]]; then
   pointing_za_deg=$5
fi

re_generate_model=1
if [[ -n "$6" && "$6" != "-" ]]; then      
   re_generate_model=$6
fi

root_opt="-q -l"
if [[ -n "$7" && "$7" != "-" ]]; then
   root_opt=$7
fi

model_options=""
if [[ -n "$9" && "$9" != "-" ]]; then
   model_options=$9
fi

start_offset=0
if [[ -n "${10}" && "${10}" != "-" ]]; then   
   start_offset=${10}
fi

beam_type=0 # 0-eda, 1-mwa
if [[ -n "${11}" && "${11}" != "-" ]]; then
   beam_type=${11}
fi

lst_start=13 # default was 15 - 20 range 
if [[ -n "${12}" && "${12}" != "-" ]]; then
   lst_start=${12}
fi

lst_end=17 # default was 15 - 20 range 
if [[ -n "${13}" && "${13}" != "-" ]]; then
   lst_end=${13}
fi

fit_time_offset=0 # was 1 for the EDA fits
if [[ -n "${14}" && "${14}" != "-" ]]; then
   fit_time_offset=${14}
fi

model_dir=""
if [[ -n "${15}" && "${15}" != "-" ]]; then
   model_dir=${15}
fi

feko_file=
if [[ -n "${16}" && "${16}" != "-" ]]; then
   feko_file=${16}
fi


# get freq. :
freq_start=`echo ${total_power_file}  | cut -b 13-15`
freq_end=`echo ${total_power_file}  | cut -b 17-19`
freq_center=`echo "$freq_start $freq_end" | awk '{print ($1+$2)/2.00;}'`
echo "Freq. range = $freq_start - $freq_end MHz -> freq center = $freq_center"

# get start uxtime :
ux_start=`head -2 ${total_power_file} | tail -1 | awk '{printf("%d\n",strtonum($1));}'`
ux_start=$(($ux_start + $start_offset))
echo "ux_start = $ux_start"
ux_end=`tail -1 ${total_power_file} | awk '{printf("%d\n",strtonum($1));}'`
echo "ux_end = $ux_end"
interval=$(($ux_end-$ux_start))

if [[ $re_generate_model -gt 0 ]]; then
   mkdir -p ${outdir}
   cd ${outdir}
   
   echo "simulate_total_power_vs_time! ${freq_center} ${ux_start} ${step} ${pointing_az_deg} ${pointing_za_deg} ${skymodel} ${freq_start} ${freq_end} \"${model_options}\" $interval - ${beam_type} ${feko_file}"
   simulate_total_power_vs_time! ${freq_center} ${ux_start} ${step} ${pointing_az_deg} ${pointing_za_deg} ${skymodel} ${freq_start} ${freq_end} "${model_options}" $interval - ${beam_type} ${feko_file}
   
   if [[ -n ${model_dir} && ${model_dir} != "-" && ${model_dir} != ${outdir} ]]; then
      dir_final="${model_dir}MODEL"
      
      echo "mv ${outdir} ${model_dir}/${dir_final}"
      mv ${outdir} ${model_dir}/${dir_final}
      
      echo "ln -s ${model_dir}/${dir_final} ${outdir}"
      ln -s ${model_dir}/${dir_final} ${outdir} 
   fi
      
   cd -
else
   echo "WARNING : re-generating of the model is not required"   
fi
# model_file=`ls ${outdir}/tant_vs_lst_*.out | tail -1`
model_file=`ls ${outdir}/tant_vs_lst_*.out | grep -v sigma | tail -1`

if [[ -s ${model_file} ]]; then
   echo "Model file $model_file exists OK"
else
   echo "ERROR : model file $model_file does not exist -> cannot continue"
   exit;
fi

total_power_file_vs_lst=${total_power_file%%.txt}_vs_LST.txt
total_power_file_vs_lst_tmp=${total_power_file%%.txt}_vs_LST.tmp
# lst_start=`ux2sid ${ux_start} mro | awk '{print $8;}'`
# echo "lst_start = $lst_start"
# awk -v lst_start=${lst_start} -v ux_start=-1 '{if($1!="#"){ux=$1;if(ux_start<0){ux_start=$1;}diff=(ux-ux_start);lst=lst_start+(diff/3600.00);while(lst>24){lst=lst-24;}print lst" "$2" "ux;}}' ${total_power_file} | sort -n > ${total_power_file_vs_lst_tmp}
# ux2lst! ${total_power_file} > ${total_power_file_vs_lst}
ux2sid_file ${total_power_file} ${total_power_file_vs_lst_tmp} mro
cat ${total_power_file_vs_lst_tmp} | sort -n > ${total_power_file_vs_lst}

model_file=`ls ${outdir}/tant_vs_lst_*.out | grep -v sigma | tail -1`
norm_rawdata_to_model! ${model_file} ${total_power_file_vs_lst}

echo "Plotting data file : ${total_power_file_vs_lst}"
echo "Model file         : ${model_file}"
# root ${root_opt} "fit_trcv_from_data_vs_time.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\")"

# use super script with Time offset fitting too :
# _ALL removed 
mkdir -p images/${freq_center}.00/
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",${lst_start},${lst_end},${fit_time_offset},\"trcv_vs_freq_LST_${lst_start}_${lst_end}_hours.txt\")"

# awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_${lst_start}_${lst_end}_hours.txt >  trcv_vs_freq_LST_${lst_start}_${lst_end}_hours_ERR.txt
# root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_${lst_start}_${lst_end}_hours_ERR.txt\")"

# 15-20 too :
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",15,20,${fit_time_offset},\"trcv_vs_freq_LST_15_20_hours.txt\")"
# awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_15_20_hours.txt >  trcv_vs_freq_LST_15_20_hours_ERR.txt
# root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_15_20_hours_ERR.txt\")"


# 11-14 - might be optimally linear for the MWA data :
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",11,14,${fit_time_offset},\"trcv_vs_freq_LST_11_14_hours.txt\")"

# 12-14 - might be optimally linear for the MWA data :
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",12,14,${fit_time_offset},\"trcv_vs_freq_LST_12_14_hours.txt\")"
# awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_15_20_hours.txt >  trcv_vs_freq_LST_15_20_hours_ERR.txt
# root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_15_20_hours_ERR.txt\")"

# 14-17 - might be optimally linear for the MWA data :
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",14,17,${fit_time_offset},\"trcv_vs_freq_LST_14_17_hours.txt\")"

# 20-22 - might be optimally linear for the MWA data :
root ${root_opt} "fit_trcv_from_data_vs_time_NEW_binned_offset.C+(\"${total_power_file_vs_lst}\",\"${model_file}\",\"power\",${freq_center},\"images/${freq_center}.00/\",20,22,${fit_time_offset},\"trcv_vs_freq_LST_20_22_hours.txt\")"


# check LST vs LOCALTIME :
ux_start=`head -2 ${total_power_file}  | tail -1 | awk '{if($1!="#"){printf("%d",$1);}}'`
plot_lst_vs_local! $ux_start - "${root_opt}"


# plot data vs model
max_diff_hours=0.016
data_vs_model_file=model_vs_data_${freq_center}MHz.out 
echo "calctxtfile ${model_file} ${total_power_file_vs_lst} ${data_vs_model_file}  -o 'v' -x ${max_diff_hours} -q use_once=1"
calctxtfile ${model_file} ${total_power_file_vs_lst} ${data_vs_model_file}  -o 'v' -x ${max_diff_hours} -q use_once=1
echo "root images/${freq_center}.00"
root ${root_opt} "plot_power_vs_model_time.C(\"${data_vs_model_file}\",2,0,NULL,\"images/${freq_center}.00/\")"

# calctxtfile total_power_080_090_MHz_MODEL/tant_vs_uxtime_1480723755_300sec_85MHz.out total_power_080_090_MHz.txt test_unixtime.out -o 'v'
data_vs_model_file_ux=model_vs_data_${freq_center}MHz_uxtime.out 
model_file_vs_uxtime=`ls ${outdir}/tant_vs_uxtime_*.out | grep -v sigma | tail -1`
echo "calctxtfile ${model_file_vs_uxtime} ${total_power_file} ${data_vs_model_file_ux} -o 'v' -x ${max_diff_hours} -q use_once=1"
calctxtfile ${model_file_vs_uxtime} ${total_power_file} ${data_vs_model_file_ux} -o 'v' -x ${max_diff_hours} -q use_once=1
echo "root images/${freq_center}.00"
root ${root_opt} "plot_power_vs_model_time.C(\"${data_vs_model_file_ux}\",2,1,NULL,\"images/${freq_center}.00/\")"

# TEST:
echo "plot_all_trcv_from_power_vs_model.sh model_vs_data_${freq_center}MHz - \"${root_opt}\""
plot_all_trcv_from_power_vs_model.sh model_vs_data_${freq_center}MHz - "${root_opt}"


