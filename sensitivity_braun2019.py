#!/usr/bin/python
from __future__ import print_function
import sys
import math
import numpy
from scipy.interpolate import interp1d




# Based on : /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{print $4" "$1;}' lfaa_requirements_20170330.rw_email
# /home/msok/Desktop/EDA2/papers/2020/M_Sokolowski/EuCAP/M_Sokolowski/logbook/20200929_lfaa_specifications.odt
def lfaa_ska( freq_mhz, interpolation_kind='cubic' ) :
   freq_mhz_array = [ 50, 55,  80, 110, 140, 160, 220, 280, 340, 345, 350 ] # /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{printf("%s, ",$4);}' lfaa_requirements_20170330.rw_email
   station_aot    = [ 68, 70, 232, 531, 588, 610, 614, 576, 522, 515, 516 ] # /home/msok/Desktop/EDA/data/sensitivity/20161207_eda_to_mwa/3c444_evening/full_band_sensitivity$ awk '{printf("%s, ",$1);}' lfaa_requirements_20170330.rw_email
   
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
   return aot_ret
   

def lfaa_sefd( freq_mhz, n_stations=512, interpolation_kind='cubic' ) :
   aot = lfaa_ska( freq_mhz )
   aot_one = aot / n_stations
   
   return ((2.00*1380.00)/aot_one)
   
def aot2sefd( aot ) :
   return ((2.00*1380.00)/aot)   
   
   
def lfaa_per_station( freq_mhz , n_stations=512, interpolation_kind='cubic' ) :
   aot_ska = lfaa_ska( freq_mhz, interpolation_kind=interpolation_kind )
   
   return ( aot_ska / n_stations )


if __name__ == "__main__":
    freq_mhz = 160.00
    if len(sys.argv) >= 1 :
       freq_mhz = float( sys.argv[1] )

    interpolation_kind='cubic'
    if len(sys.argv) >= 2 :
       interpolation_kind = sys.argv[2]

    aot_ska = lfaa_ska( freq_mhz , interpolation_kind=interpolation_kind )
    sefd_ska = aot2sefd( aot_ska )
    
    aot_ska_station = lfaa_per_station( freq_mhz , interpolation_kind=interpolation_kind )
    sefd_ska_station = aot2sefd( aot_ska_station )
    
    print("A/T_ska         ( %.2f MHz) = %.4f m^2/K  or SEFD_ska     = %.2f Jy" % (freq_mhz,aot_ska,sefd_ska))
    print("A/T_ska_station ( %.2f MHz) = %.4f m^2/K  or SEFD_station = %.2f Jy" % (freq_mhz,aot_ska_station,sefd_ska_station))
       
