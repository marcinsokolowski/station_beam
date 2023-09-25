#!/bin/bash

basename=model_vs_data_85MHz
if [[ -n "$1" && "$1" != "-" ]]; then
   basename=$1
fi
freq=`echo $basename | awk '{freq_last=substr($1,17,1);len=3;if(freq_last=="M"){len=2;}freq=substr($1,15,len);print freq;}'`
echo "Freq = $freq -> subdir = images/${freq}.00"

outfile=${basename}_FIT.out
if [[ -n "$2" && "$2" != "-" ]]; then
   outfile=$2
fi
rm -f ${outfile}

root_opt="-l -q"
if [[ -n "$3" && "$3" != "-" ]]; then
   root_opt=$3
fi

mkdir -p images/
split_by_lst! ${basename}.out
ls ${basename}_[0-9]*_[0-9]*.out| grep -v sigma > ${basename}.list
root ${root_opt} "plot_power_vs_model.C(\"${basename}.list\",\"images/${freq}.00/\")"

for infile in `ls ${basename}_[0-9]*_[0-9]*.out | grep -v sigma`
do
   if [[ -s ${infile} ]]; then
      root ${root_opt} "fit_trcv_from_data_vs_model.C(\"${infile}\",\"line\",\"${outfile}\")"
   else
      echo "ERROR : could not execute fit_trcv_from_data_vs_model.C(\"${infile}\",\"line\",\"${outfile}\")"
      echo "ERROR : input file ${infile} does not exist"
   fi
done

if [[ -s ${basename}.out ]]; then
   root ${root_opt} "fit_trcv_from_data_vs_model.C(\"${basename}.out\",\"line\",\"${outfile}\")"
else
   echo "ERROR : could not execute fit_trcv_from_data_vs_model.C(\"${basename}.out\",\"line\",\"${outfile}\")"
   echo "ERROR : input file ${basename}.out does not exist"
fi

awk '{print ($2+$3)/2.00" "$4;}' ${outfile} > ${basename}_Trcv_vs_LST.out
awk '{print ($2+$3)/2.00" "$6;}' ${outfile} > ${basename}_Gain_vs_LST.out

if [[ -s ${basename}_Gain_vs_LST.out ]]; then
   root ${root_opt} "plot_gain_vs_lst.C(\"${basename}_Gain_vs_LST.out\",\"images/${freq}.00/\")"
else
   echo "ERROR : could not execute plot_gain_vs_lst.C(\"${basename}_Gain_vs_LST.out\",\"images/${freq}.00/\")"
   echo "ERROR : input file ${basename}_Gain_vs_LST.out does not exist"
fi

if [[ -s ${basename}_Trcv_vs_LST.out ]]; then
   root ${root_opt} "plot_trcv_vs_lst.C(\"${basename}_Trcv_vs_LST.out\",\"images/${freq}.00/\")"
else
   echo "ERROR : could not execute plot_trcv_vs_lst.C(\"${basename}_Trcv_vs_LST.out\",\"images/${freq}.00/\")"
   echo "ERORR : input file ${basename}_Trcv_vs_LST.out does not exist"
fi

