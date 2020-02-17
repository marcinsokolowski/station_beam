#!/bin/bash

creator="msok"
code_version="version0.00"
station_type=0 # EDA2

for txtfile in `ls *.txt`
do
  sqlfile=${txtfile%%txt}sql
  b=${txtfile%%.txt}
  polarisation=`echo $b | awk '{last2=length($1);print substr($1,last2);}'`
  
  header_line=`grep HEADER $txtfile | tail -1`
  header_count=`grep HEADER $txtfile | wc -l`
  
  gpstime=`echo $header_line | awk '{print int($6);}'`
  ux=`gps2ux! $gpstime`
  lst=`ux2sid $ux | awk '{print $8;}'`
  az=`echo $header_line | awk '{print $10;}'`
  za=`echo $header_line | awk '{print $15;}'`
  
  echo 
  echo "$header_line / $header_count -> $sqlfile"
  echo "   gps = $gpstime, az = $az, za = $za -> ux=$ux , lst=$lst , polarisation = $polarisation"
  

# FREQ[MHz] Sensitivity[m^2/K] T_sys[K] A_eff[m^2] T_rcv[K] T_ant[K] M=(1-|G_a|^2)xF
# 49.92000000 0.08786082 11925.35 1047.77075195 3851.41133725 8073.93501201 1.00000000

  awk -v pol=${polarisation} -v header_count=${header_count} -v cnt=0 -v start=0 -v az=${az} -v za=${za} -v lst=${lst} -v ux=${ux} -v gps=${gpstime} -v station_type=${station_type} -v creator=${creator} -v code_version=${code_version} \
  '{
     if( $2 == "HEADER" ){
        cnt += 1;

        if( cnt == header_count ){
           start=1;
        }
     }else{
        if( $1 != "#" ){
           if( start > 0 ){
               printf("INSERT INTO Sensitivity (polarisation,azim_deg,za_deg,frequency_mhz,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,creator,code_version) VALUES (\"%s\"%.4f,%.4f,%.2f,%.8f,%.2f,%.2f,%.8f,%.2f,%.4f,%.2f,%.2f,%d,\"%s\",\"%s\");\n",pol,az,za,$1,lst,ux,gps,$2,$3,$4,$5,$6,station_type,creator,code_version);
           }
         }
     }
  }'  $txtfile > $sqlfile

done
