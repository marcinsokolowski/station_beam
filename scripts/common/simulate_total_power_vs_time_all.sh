#!/bin/bash

step=300
if [[ -n "$1" && "$1" != "-" ]]; then
   step=$1
fi         

pointing_az_deg=0
if [[ -n "$2" && "$2" != "-" ]]; then
   pointing_az_deg=$2
fi
          
pointing_za_deg=0
if [[ -n "$3" && "$3" != "-" ]]; then
   pointing_za_deg=$3
fi
              
re_generate_model=1
if [[ -n "$4" && "$4" != "-" ]]; then
   re_generate_model=$4
fi

root_opt="-q -l"
if [[ -n "$5" && "$5" != "-" ]]; then
   root_opt=$5
fi

skymodel=MODEL
if [[ -n "$6" && "$6" != "-" ]]; then
   skymodel=$6
fi

list=`ls total_power_[0-9]*_[0-9]*_MHz.txt`
if [[ -n "$7" && "$7" != "-" ]]; then
   list=$7
fi
   
# if [[ -n "$8" && "$8" != "-" ]]; then

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

# outdir=${total_power_file%%.txt}_${skymodel}/
model_dir=""
if [[ -n "${14}" && "${14}" != "-" ]]; then
   model_dir=${14}
   # re_generate_model=0
   # echo "model_dir set to $model_dir -> re_generate_model -> 0"
fi

feko_file=
if [[ -n "${15}" && "${15}" != "-" ]]; then
   feko_file=${15}
fi

echo "#######################################################################"
echo "PARAMETERS:"
echo "#######################################################################"
echo "List of files     = $list"
echo "beam type         = $beam_type"
echo "LST-range of fit  = $lst_start - $lst_end [hours]"
echo "Model dir         = $model_dir"
echo "re_generate_model = $re_generate_model"
echo "feko_file         = $feko_file"
echo "#######################################################################"
            
# backup of older file
dtm=`date +%Y%m%d`
mv trcv_vs_freq.txt trcv_vs_freq.txt.${dtm}

for total_power_file in `echo $list`
do
   re_generate_model_freq=$re_generate_model   
   model_dir_freq="-"
   
   if [[ -n "${model_dir}" && "${model_dir}" != "-" ]]; then
      model_dir_freq=${total_power_file%%.txt}_${skymodel}
      echo "model_dir = $model_dir -> model_dir_freq = $model_dir_freq , re_generate_model_freq = $re_generate_model_freq ( re_generate_model = $re_generate_model )" 
      
      model_sub_dir=`echo $total_power_file | cut -b 1-23 | awk -v model=${skymodel} '{print $1"_"model}'`
      full_model_dir=${model_dir}/${model_sub_dir}
      if [[ -d ${full_model_dir} ]]; then
         echo "ln -s ${full_model_dir} ${model_dir_freq}"
         ln -s ${full_model_dir} ${model_dir_freq}
         re_generate_model_freq=0
      else
         echo "WARNING : model directory ${full_model_dir} does not exist -> re_generate_model_freq = $re_generate_model_freq ( re_generate_model = $re_generate_model )"         
      fi
      
   else
      echo "No model dir provided and re_generate_model_freq = $re_generate_model_freq ( re_generate_model = $re_generate_model ) -> appropriate action will be taken"
   fi

   echo "simulate_total_power_vs_time.sh ${total_power_file} ${model_dir_freq} ${step} ${pointing_az_deg} ${pointing_za_deg} ${re_generate_model_freq} \"${root_opt}\" ${skymodel} \"${model_options}\" - ${beam_type} $lst_start $lst_end - \"${model_dir}\" ${feko_file}"
   simulate_total_power_vs_time.sh ${total_power_file} ${model_dir_freq} ${step} ${pointing_az_deg} ${pointing_za_deg} ${re_generate_model_freq} "${root_opt}" ${skymodel} "${model_options}" - ${beam_type} $lst_start $lst_end - "${model_dir}" ${feko_file}
done

# plot selected LST range :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_${lst_start}_${lst_end}_hours.txt >  trcv_vs_freq_LST_${lst_start}_${lst_end}_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_${lst_start}_${lst_end}_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_${lst_start}_${lst_end}_hours_ERR.txt\",1)"


# TODO : verify if these two things work ok :
awk '{if($1!="#" && (strtonum($1)<130 || strtonum($1)>140)){print $0}}' trcv_vs_freq.txt > trcv_vs_freq.selected
root -l -q "plot_tlna_MHz.C(\"trcv_vs_freq.selected\",NULL,0,2000)"


tail --lines=1 model_vs_data_*MHz_FIT.out | grep -v FIT | sort -n  | awk '{if(length($0)>0){print $1" "$4;}}' > trcv_vs_freq.data_vs_model_all
root -l -q "plot_tlna_MHz.C(\"trcv_vs_freq.data_vs_model_all\",NULL,0,2000)"

echo "cp trcv_vs_freq.txt trcv_vs_freq_LST_${lst_start}_${lst_end}.txt"
cp trcv_vs_freq.txt trcv_vs_freq_LST_${lst_start}_${lst_end}.txt

echo "cp trcv_vs_freq.selected trcv_vs_freq_LST_${lst_start}_${lst_end}.selected"
cp trcv_vs_freq.selected trcv_vs_freq_LST_${lst_start}_${lst_end}.selected


echo "Plotting extra LST ranges :"
date

# Additional plots using different LST ranges :
# plot 15-20 LST range too :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_15_20_hours.txt >  trcv_vs_freq_LST_15_20_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_15_20_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_15_20_hours_ERR.txt\",1)"

# plot 12-14 LST range too :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_12_14_hours.txt >  trcv_vs_freq_LST_12_14_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_12_14_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_12_14_hours_ERR.txt\",1)"

# plot 11-14 LST range too :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_11_14_hours.txt >  trcv_vs_freq_LST_11_14_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_11_14_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_11_14_hours_ERR.txt\",1)"

# plot 14-17 LST range too :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_14_17_hours.txt >  trcv_vs_freq_LST_14_17_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_14_17_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_14_17_hours_ERR.txt\",1)"

# plot 20-22 LST range too :
awk '{print $1" 5 "$2" "$3;}' trcv_vs_freq_LST_20_22_hours.txt >  trcv_vs_freq_LST_20_22_hours_ERR.txt
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_20_22_hours_ERR.txt\")"
root ${root_opt} "plot_trcv_vs_freq_with_err.C(\"trcv_vs_freq_LST_20_22_hours_ERR.txt\",1)"


