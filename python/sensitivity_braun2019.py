#!/usr/bin/python
from __future__ import print_function
import sys
import math
import numpy
from scipy.interpolate import interp1d




# Table 9, column 2 in Braun (2019) : braun_2019.txt
# awk '{printf("%d, ",$1);}' braun_2019.txt
# awk '{printf("%d, ",$2);}' braun_2019.txt
def skalow_braun2019( freq_mhz, interpolation_kind='cubic' ) :
   freq_mhz_array = [ 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350 ] # awk '{printf("%d, ",$1);}' braun_2019.txt
   station_aot    = [ 0.102, 0.255, 0.420, 0.580, 0.781, 0.862, 0.885, 0.924, 0.965, 1.052, 1.119, 1.159, 1.171, 1.190, 1.195, 1.223, 1.235, 1.263, 1.226, 1.234, 1.282, 1.304, 1.297, 1.306, 1.314, 1.289, 1.294, 1.327, 1.286, 1.204, 1.156 ] # awk '{printf("%.3f, ",$2);}' braun_2019.txt
   
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
   
   aot_ret = aot_func( freq_mhz )*512 # 512 stations 
   return aot_ret
   

def lfaa_sefd( freq_mhz, n_stations=512, interpolation_kind='cubic' ) :
   aot = skalow_braun2019( freq_mhz )
   aot_one = aot / n_stations
   
   return ((2.00*1380.00)/aot_one)
   
def aot2sefd( aot ) :
   return ((2.00*1380.00)/aot)   
   
   
def lfaa_per_station( freq_mhz , n_stations=512, interpolation_kind='cubic' ) :
   aot_ska = skalow_braun2019( freq_mhz, interpolation_kind=interpolation_kind )
   
   return ( aot_ska / n_stations )


if __name__ == "__main__":
    freq_mhz = 160.00
    if len(sys.argv) >= 1 :
       freq_mhz = float( sys.argv[1] )

    interpolation_kind='cubic'
#    if len(sys.argv) >= 3 :
#       interpolation_kind = sys.argv[2]

    aot_ska = skalow_braun2019( freq_mhz , interpolation_kind=interpolation_kind )
    sefd_ska = aot2sefd( aot_ska )
    
    aot_ska_station = lfaa_per_station( freq_mhz , interpolation_kind=interpolation_kind )
    sefd_ska_station = aot2sefd( aot_ska_station )
    
    print("A/T_ska         ( %.2f MHz) = %.4f m^2/K  or SEFD_ska     = %.2f Jy" % (freq_mhz,aot_ska,sefd_ska))
    print("A/T_ska_station ( %.2f MHz) = %.4f m^2/K  or SEFD_station = %.2f Jy" % (freq_mhz,aot_ska_station,sefd_ska_station))
       
