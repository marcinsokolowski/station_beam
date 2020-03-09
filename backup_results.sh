#!/bin/bash

do_rm=0
if [[ -n "$1" && "$1" != "-" ]]; then
   do_rm=$1
fi

dtm=`date +%Y%m%d_%H%M%S`
mkdir -p backup/${dtm}/

echo "cp 1*.txt backup/${dtm}/"
cp 1*.txt backup/${dtm}/

echo "mv out backup/${dtm}/out.${dtm}"
mv out backup/${dtm}/out.${dtm}

cd backup/${dtm}/
echo "nohup gzip out.${dtm} > gzip.out 2>&1 &"
nohup gzip out.${dtm} > gzip.out 2>&1 &
cd -


ls -ltr 1*.txt
ls -tr 1*.txt | tail --lines=2 | awk 'BEGIN{rm_line="rm -f ";}{rm_line = rm_line " " $1;}END{print rm_line;}'


if [[ $do_rm -gt 0 ]]; then
   rm_cmd=`ls -tr 1*.txt | tail --lines=2 | awk 'BEGIN{rm_line="rm -f ";}{rm_line = rm_line " " $1;}END{print rm_line;}'`
   
   echo "Executing : $rm_cmd"
   $rm_cmd
else
   echo "WARNING : automatic removing of last file is not required"
fi
