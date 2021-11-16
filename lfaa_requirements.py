#!/usr/bin/python
from __future__ import print_function
import sys
import math
import numpy
from scipy.interpolate import interp1d


# NEW (v12) : /home/msok/Desktop/SKA/papers/2020/M_Sokolowski/Sensitivity_DB/SUBMISSION/REVIEW/references/Sensitivity_v12_or_whatever.txt
#             ( every 10 MHz and includes average above elevation 45 deg)
#             cd /home/msok/Desktop/SKA/papers/2020/M_Sokolowski/Sensitivity_DB/SUBMISSION/REVIEW/references
#             awk '{if($1!="#"){printf("%.3f, ",$1);}}' Sensitivity_v12_or_whatever.txt
#             awk '{if($1!="#"){printf("%.3f, ",$2);}}' Sensitivity_v12_or_whatever.txt
#
#
# OLD : /home/msok/Desktop/EDA2/papers/2020/M_Sokolowski/EuCAP/M_Sokolowski/logbook/20200929_lfaa_specifications.odt
#          (only at zenith and less frequencies)
#        Based on : /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{print $4" "$1;}' lfaa_requirements_20170330.rw_email
#        /home/msok/Desktop/EDA2/papers/2020/M_Sokolowski/EuCAP/M_Sokolowski/logbook/20200929_lfaa_specifications.odt
def lfaa_ska( freq_mhz, interpolation_kind='cubic' , version="v12" , subversion="AvgAboveElev45deg" ) : # subversion : AvgAboveElev45deg or Zenith
   # OLD was a default 
   freq_mhz_array = [ 50, 55,  80, 110, 140, 160, 220, 280, 340, 345, 350 ] # /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{printf("%s, ",$4);}' lfaa_requirements_20170330.rw_email
   station_aot    = [ 68, 70, 232, 531, 588, 610, 614, 576, 522, 515, 516 ] # /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{printf("%s, ",$1);}' lfaa_requirements_20170330.rw_email
   
   if version == "v12" or version == "12" :
      if subversion.lower() == "avgaboveelev45deg" :
         freq_mhz_array = [ 50.000, 60.000,  70.000,  80.000,  90.000, 100.000, 110.000, 120.000, 130.000, 140.000, 150.000, 160.000, 170.000, 180.000, 190.000, 200.000, 210.000, 220.000, 230.000, 240.000, 250.000, 260.000, 270.000, 280.000, 290.000, 300.000, 310.000, 320.000, 330.000, 340.000, 350.000 ]
         station_aot    = [ 56.000, 68.000, 115.000, 190.000, 280.000, 368.000, 435.000, 469.000, 479.000, 482.000, 490.000, 500.000, 508.000, 511.000, 512.000, 511.000, 508.000, 504.000, 499.000, 493.000, 488.000, 483.000, 477.000, 472.000, 468.000, 463.000, 457.000, 450.000, 441.000, 428.000, 423.000 ]
      else: # default is zenith :
         freq_mhz_array = [ 50.000, 60.000,  70.000,  80.000,  90.000, 100.000, 110.000, 120.000, 130.000, 140.000, 150.000, 160.000, 170.000, 180.000, 190.000, 200.000, 210.000, 220.000, 230.000, 240.000, 250.000, 260.000, 270.000, 280.000, 290.000, 300.000, 310.000, 320.000, 330.000, 340.000, 350.000 ]
         station_aot    = [ 68.000, 83.000, 141.000, 232.000, 342.000, 449.000, 531.000, 572.000, 584.000, 588.000, 598.000, 610.000, 619.000, 624.000, 625.000, 623.000, 619.000, 614.000, 608.000, 602.000, 595.000, 588.000, 582.000, 576.000, 570.000, 565.000, 558.000, 549.000, 537.000, 522.000, 516.000 ]         
         
      print("DEBUG : used sensitivity requirements %s/%s" % (version,subversion))
         
   
   aot_func = None
   
   last = len(freq_mhz_array)-1 
   if freq_mhz >= freq_mhz_array[0] and freq_mhz <= freq_mhz_array[last] :         
      aot_func = interp1d( freq_mhz_array , station_aot , kind=interpolation_kind )
   else :
      print("WARNING : freq_mhz outside the SKA Specs region %.2f - %.2f MHz" % (freq_mhz_array[0],freq_mhz_array[last]))
      if freq_mhz < freq_mhz_array[0] :
         return station_aot[0]
      else :
         return station_aot[last]

   if aot_func is None :
      print("ERROR : aot_func is None -> no A/T value returned")
      return None         
   
   aot_ret = aot_func( freq_mhz )
   print("DEBUG : lfaa_ska( %.2f MHz) returning %.4f [m^2/K]" % (freq_mhz,aot_ret))
   return aot_ret
   

def lfaa_sefd( freq_mhz, n_stations=512, interpolation_kind='cubic' ) :
   aot = lfaa_ska( freq_mhz )
   aot_one = aot / n_stations
   
   return ((2.00*1380.00)/aot_one)
   
def aot2sefd( aot ) :
   return ((2.00*1380.00)/aot)   
   
   
def lfaa_per_station( freq_mhz , n_stations=512, interpolation_kind='cubic', version="v12" , subversion="AvgAboveElev45deg" ) :
   aot_ska = lfaa_ska( freq_mhz, interpolation_kind=interpolation_kind, version=version, subversion=subversion )
   
   return ( aot_ska / n_stations )


if __name__ == "__main__":
#    print("DEBUG : number of parameters = %d , they are : %s %s" % (len(sys.argv),sys.argv[0],sys.argv[1]))

    freq_mhz = 160.00
    if len(sys.argv) >= 1 :
       freq_mhz = float( sys.argv[1] )

    interpolation_kind='cubic'
    if len(sys.argv) >= 3 :
       interpolation_kind = sys.argv[2]
       
    version="v12"
    if len(sys.argv) >= 4 :   
       version = sys.argv[3]
       
    subversion="AvgAboveElev45deg"   
    if len(sys.argv) >= 5 :
       subversion = sys.argv[4]
    
    print("#######################################################################")   
    print("PARAMETERS :")
    print("#######################################################################")
    print("freq_mhz = %.4f [MHz]" % (freq_mhz))
    print("interpolation_kind = %s" % (interpolation_kind))
    print("version  = %s/%s" % (version,subversion))
    print("#######################################################################")
       

    aot_ska = lfaa_ska( freq_mhz , interpolation_kind=interpolation_kind, version=version, subversion=subversion )
    sefd_ska = aot2sefd( aot_ska )
    
    aot_ska_station = lfaa_per_station( freq_mhz , interpolation_kind=interpolation_kind, version=version, subversion=subversion )
    sefd_ska_station = aot2sefd( aot_ska_station )
    
    print("A/T_ska         ( %.2f MHz) = %.4f m^2/K  or SEFD_ska     = %.2f Jy" % (freq_mhz,aot_ska,sefd_ska))
    print("A/T_ska_station ( %.2f MHz) = %.4f m^2/K  or SEFD_station = %.2f Jy" % (freq_mhz,aot_ska_station,sefd_ska_station))
       
