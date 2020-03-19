#!/usr/bin/python
import sys
import math
import numpy
from scipy.interpolate import interp1d

# see ~/Desktop/EDA/loogbook/Cath/Budi_trcv and ~/Desktop/EDA/loogbook/Cath/eda_aeff_cath.odt
# Budi's e-mail on 2016-12-20 
def trcv_budi( freq_mhz ):
   t_lna =  +724.02940227341059653554111719131469726562500000000000 -18.91813674541278089691331842914223670959472656250000*freq_mhz \
             +0.23407541076266874524591798945039045065641403198242*freq_mhz*freq_mhz -0.00233671195682318014533174199698351003462448716164*freq_mhz*freq_mhz*freq_mhz \
             +0.00002254381274811321963493990216953477556671714410*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000016409668824170788921994429236661527937712890*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000076016725276938624039072744188781996510507*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz - 0.00000000000210443395476647521171394324012829569895*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000317128431524229797923526252249329818*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz - 0.00000000000000000179318498778901811500861463716146*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000034683893039093007945200578004*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             - 0.00000000000000000000000036757731954965862822568812*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000262855743668854032145143 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000001701169748990267910934 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000001832975188922028666 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000010411258968305159 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000036864325213161 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000000061278189979 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000000000147452601 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000000000000000000001736900 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   print "INFO : using Budi's data for T_rcv ( %.2f MHz) = %.2f" % (freq_mhz,t_lna)

   return t_lna   

# fitted from POWER vs. MODEL data in LST range 18-20 hours only see Ill. Illustration 6 in eda_drift_scan_data_201612.odt for details:
def trcv_fit_data_vs_mode( freq_mhz ):
   t_lna = +40641.44874588610400678589940071105957031250000000000000  -1151.44119960312536932178772985935211181640625000000000*freq_mhz \
           +10.28736969420266333941071934532374143600463867187500 * freq_mhz*freq_mhz -0.00920101512017338675486488597243805998004972934723 * freq_mhz*freq_mhz*freq_mhz -0.00024102523369100168010650819816476086998591199517 * freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           -0.00000019651325627217189202790260330738769667391352 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000820101610522132312859120128458942367011275 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           +0.00000000001047686030411538916735193183367932131915 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000000024494646525009589774756356651373232348 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           +0.00000000000000052438247226051691783257993415371231 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   return t_lna

# Fitted - see :
# cd /home/msok/Desktop/EDA/data/2016-07/20160703
# root [4] .x plot_tlna_MHz.C("trcv_vs_freq.selected",NULL,0,5000)
# t_lna =  +36853.38919790466752601787447929382324218750000000000000 * 1 -1377.64458069924739902489818632602691650390625000000000 * 1*freq_mhz +20.37196908756931890138730523176491260528564453125000 * 1*freq_mhz*freq_mhz -0.13515251806813355361924777753301896154880523681641 * 1*freq_mhz*freq_mhz*freq_mhz +0.00022215700342832035230414278181854115246096625924 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000169842384409759492976880849990362065682347747 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000464297999395697167928294570157440879221156 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000003336389317067374617426693716906293402102 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000000017598078793649799250229898970650520328 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000000000023477806740472191603243644100388291 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz;
# Cubic interpolation of data points see  T_rcv from lightcurve fits see Ill.1 and Ill.2 in trcv_from_20160703_drift_scan.odt for details
def trcv_fit_lightcurve_20160703_cubic( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' trcv_vs_freq.selected
   # awk '{printf("%.2f,",$2);}' trcv_vs_freq.selected
   x=[55.00,65.00,75.00,85.00,95.00,105.00,115.00,145.00,155.00,165.00,175.00,185.00,195.00,205.00,215.00,225.00,235.00]
   y=[3831.09,1744.61,1022.83,695.13,489.77,340.59,232.78,83.88,59.21,45.78,35.69,30.60,27.26,30.52,27.69,28.39,93.70]
      
   tlna_cubic = interp1d(x, y, kind='cubic')
   
   return tlna_cubic(freq_mhz)

# EDA/loogbook/haslam_vs_angelica.odt Ill. 24 red curve
# /home/msok/Desktop/EDA/data/2016-07/20160703/ANGELICA/trcv_vs_freq.selected
def trcv_angelica_data_vs_time( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' ANGELICA/trcv_vs_freq.selected   
   # awk '{printf("%.2f,",$2);}' ANGELICA/trcv_vs_freq.selected
   x=[55.00,65.00,75.00,85.00,95.00,105.00,125.00,155.00,175.00,205.00,225.00]
   y=[3322.52,1608.02,1007.08,824.80,638.77,232.16,245.80,68.14,47.27,33.88,27.79]
   
   tlna_cubic = interp1d(x, y, kind='cubic')
   
   return tlna_cubic(freq_mhz)
   
def trcv_angelica_data_vs_time_powerlawfit( freq_mhz ) :
   t_rcv = 384.3 * math.pow((100.00/freq_mhz),3.571)
   
   if freq_mhz >= 160 and False :      
      t_rcv = trcv_budi(freq_mhz)/corr_factor_adrian(freq_mhz)
      
   return t_rcv   
      

# EDA/loogbook/haslam_vs_angelica.odt Ill. 24 green curve
# /home/msok/Desktop/EDA/data/2016-07/20160703/HASLAM/trcv_vs_freq.data_vs_model_all
def trcv_haslam_data_vs_model( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' HASLAM/trcv_vs_freq.data_vs_model_all
   # awk '{printf("%.2f,",$2);}' HASLAM/trcv_vs_freq.data_vs_model_all
   x=[55.00,65.00,75.00,85.00,95.00,105.00,115.00,145.00,155.00,165.00,175.00,185.00,195.00,205.00,215.00,225.00,235.00]
   y=[3849.12,1744.54,1013.68,682.12,473.69,327.36,222.33,81.70,57.44,45.02,35.46,30.50,27.35,31.30,28.75,29.14,45.67]
   
   tlna_cubic = interp1d(x, y, kind='cubic')
   return tlna_cubic(freq_mhz)

# as in loogbook/paper/eda_lightcurve_and_trcv.odt Ill. 13:
# /home/msok/Desktop/EDA/data/2016-07/20160703/HASLAM/WithSun_and_BeamformingErrors/BIGHORNS/images/final/paper/trcv_vs_freq_lst_13-17_hours.png
# root [1] .x plotspec_freq_errxy.C("trcv_vs_freq_lst_range_13_17_hours.txt")
# Fitted : [0]*(150/x)^[1]
def trcv_from_skymodel_with_err( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too 
   # not to go below 50K (as fitting might be tricky)
#   if trcv < 50 or freq_mhz>160 :
#   if freq_mhz > 160 :
#      return 50

   print "INFO : using trcv_from_skymodel_with_err"

   if use_cubic :
      x=[55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0]
      y=[3722.55055596,1421.83480972,655.8033454,409.69506365,356.59552554,293.4326318,200.26593869,163.91940732,157.46121101,112.19530177,51.14849704,68.56693587,28.33318976,43.72596282,34.41633474,62.76422256]   
      tlna_cubic = interp1d(x, y, kind='cubic')
      trcv = tlna_cubic(freq_mhz)
   else :
      index=3.46114
      trcv = 82.4708*math.pow( (150.00/freq_mhz), index )
   
   if trcv < 50 or freq_mhz>160 :
      trcv = 50
      
   # see /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity.odt or /home/msok/Desktop/EDA/doc/T_rcv   
#   if freq_mhz >= 240 :
#      trcv = 70 

   # see /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity.odt and /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity_FINAL.odt
   if freq_mhz > 160 :
      trcv = 30.2000 + freq_mhz*0.100125
   
   return trcv


# AAVS2 : 
# see : /home/msok/Desktop/AAVS2/logbook/20200128_receiver_temperature.odt
def trcv_aavs2( freq_mhz , use_cubic=True, za_deg=None ): # default was power law fit but I wanted to be able to use cubic fit too 
   # not to go below 50K (as fitting might be tricky)
#   if trcv < 50 or freq_mhz>160 :
#   if freq_mhz > 160 :
#      return 50


   # at the moment return the same as for EDA-1 (or EDA-2)
#   return trcv_from_skymodel_with_err( freq_mhz , use_cubic )

   trcv = 33.85 + 0.0148*freq_mhz 
   
   using_cubic = False
   if freq_mhz < 100.00 :
      x=[10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,42.0,44.0,46.0,48.0,50.0,52.0,54.0,56.0,58.0,60.0,62.0,64.0,66.0,68.0,70.0,72.0,74.0,76.0,78.0,80.0,82.0,84.0,86.0,88.0,90.0,92.0,94.0,96.0,98.0,100.0,102.0,104.0,106.0,108.0,110.0,112.0,114.0,116.0,118.0,120.0,122.0,124.0,126.0,128.0,130.0,132.0,134.0,136.0,138.0,140.0,142.0,144.0,146.0,148.0,150.0,152.0,154.0,156.0,158.0,160.0,162.0,164.0,166.0,168.0,170.0,172.0,174.0,176.0,178.0,180.0,182.0,184.0,186.0,188.0,190.0,192.0,194.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,292.0,294.0,296.0,298.0,300.0,302.0,304.0,306.0,308.0,310.0,312.0,314.0,316.0,318.0,320.0,322.0,324.0,326.0,328.0,330.0,332.0,334.0,336.0,338.0,340.0,342.0,344.0,346.0,348.0,350.0,352.0,354.0,356.0,358.0,360.0,362.0,364.0,366.0,368.0,370.0,372.0,374.0,376.0,378.0,380.0,382.0,384.0,386.0,388.0,390.0,392.0,394.0,396.0,398.0,400.0,402.0,404.0,406.0,408.0,410.0]
      y=[7107100.0000,654850.0000,206830.0000,376900.0000,361870.0000,361320.0000,384580.0000,408710.0000,370030.0000,257710.0000,194070.0000,108050.0000,55808.0000,29746.0000,17932.0000,11791.0000,7751.7000,5199.5000,3493.5000,2188.8000,1356.3000,747.7300,401.0200,210.2100,118.1100,82.2620,72.2770,68.3100,64.0260,57.4240,49.0060,43.9750,40.7470,41.5500,44.0500,46.5410,48.0550,48.9760,46.7510,42.1360,39.1750,33.6010,32.7680,32.4490,32.4100,35.3490,33.8180,35.4250,36.2200,37.5340,36.2480,36.4220,35.7270,35.0100,35.8200,35.4200,37.2410,37.0950,36.3850,35.0540,33.2270,33.0210,35.3170,36.1650,38.3440,38.6990,36.5530,35.4320,34.9590,34.6170,36.0870,37.4090,38.1410,39.4630,37.1560,36.1880,36.3860,35.1380,35.4160,35.5970,37.9060,38.6590,38.3350,36.6420,37.6530,36.6690,34.2690,34.3780,33.7400,35.3170,36.9150,37.8910,37.6220,36.9020,37.7700,40.9420,35.1330,35.2780,33.5290,32.9030,34.3800,34.4300,37.6530,39.4540,39.4440,41.0410,40.7670,40.4160,40.1550,39.6900,38.7630,36.3620,35.3780,34.3560,33.9780,33.9440,33.6700,32.6710,34.1480,34.1550,34.2700,35.1280,34.9730,34.7650,34.0960,35.0130,35.2690,36.0910,36.1660,38.7710,38.3090,39.7370,39.9530,41.9870,42.0490,37.8430,40.1820,42.6590,40.5190,39.8170,40.5610,38.7800,37.9280,37.0910,37.1220,36.7460,38.1390,37.4310,37.3780,39.4030,39.4210,39.3430,40.6120,41.0690,40.2540,36.9340,40.3270,39.4150,39.3720,39.0670,38.4750,36.8540,36.4600,37.2850,37.3790,37.4560,38.3510,36.9400,38.7410,39.1510,38.8910,42.3500,41.2200,41.3130,40.9190,37.1090,39.4140,40.1410,38.1900,37.2650,37.8390,37.3610,36.5130,36.7570,37.5530,37.1360,39.0620,39.9190,42.8730,44.6690,45.2820,45.9850,47.6020,47.5390,49.5910,47.3790,49.9540,49.2800,49.2160,49.2670,47.3600]
      
      trcv_interpol = interp1d(x, y, kind='cubic')
      trcv = trcv_interpol( freq_mhz )
      using_cubic = True

   if za_deg is not None :
      print "za_deg is not None:"
      trcv_vs_za_fit_params={}
      trcv_vs_za_fit_params[80]  = numpy.array( [ 44.37610000 , 0.01785778 , 0.00073852 ] )
      trcv_vs_za_fit_params[110] = numpy.array( [ 39.43860000 , 0.02064556 , 0.00000393 ] )
      trcv_vs_za_fit_params[140] = numpy.array( [ 42.35740000 , -0.12224778 , 0.00070570 ] ) 
      trcv_vs_za_fit_params[160] = numpy.array( [ 42.44980000 , -0.06561444 , 0.00028481 ] ) 
      trcv_vs_za_fit_params[210] = numpy.array( [ 43.53470000 , -0.06032556 , 0.00094985 ] )
      trcv_vs_za_fit_params[220] = numpy.array( [ 46.33950000 , 0.07806778 , -0.00182215 ] )
      trcv_vs_za_fit_params[230] = numpy.array( [ 44.78490000 , 0.10587222 , -0.00187741 ] )
      trcv_vs_za_fit_params[280] = numpy.array( [ 50.56900000 , 0.08150778 , -0.00180215 ] )
      trcv_vs_za_fit_params[340] = numpy.array( [ 49.09320000 , -0.02795556 , 0.00065852 ] )
      trcv_vs_za_fit_params[345] = numpy.array( [ 51.47230000 , 0.01942444 , -0.00038193 ] )
      trcv_vs_za_fit_params[350] = numpy.array( [ 52.30620000 , 0.03021889 , -0.00063785 ] )
     
      if freq_mhz >= ( min(trcv_vs_za_fit_params.keys())-10 ) and freq_mhz <= ( max(trcv_vs_za_fit_params.keys())+10 ) :
         best_freq=None
         mindist=1e6
         for freq_key in trcv_vs_za_fit_params.keys() :
            dist = math.fabs( float(freq_key) - freq_mhz ) 
            if dist < mindist :
               mindist = dist
               best_freq = freq_key
         
         if best_freq is not None :
            fit_params = trcv_vs_za_fit_params[best_freq]            
            trcv_vs_za = fit_params[0] + fit_params[1]*za_deg + fit_params[2]*(za_deg**2)
            
            print "INFO : calculated trcv_aavs(%.2f MHz at za = %.2f [deg]) = %.2f [K] (vs. normal = %.2f [K]) , used fit coefficients %.4f, %.4f , %.4f" %  (freq_mhz,za_deg,trcv_vs_za,trcv,fit_params[0],fit_params[1],fit_params[2])
            trcv = trcv_vs_za                           
         else :
            print "WARNING : best_freq not found for %.2f MHz -> do not using ZA dependence" % (freq_mhz)
      else :
         print "WARNING : frequency %.2f MHz outside the range of known ZA dependence (80-350 MHz) -> do not using ZA dependence" % (freq_mhz)
   else:
      print "WARNING : za_deg is None - using T_rcv(freq) only not as a function of pointing ZA[deg]"

   print "INFO : using trcv_aavs2( %.2f MHz) = %.2f [K] (cubic = %s)" % (freq_mhz,trcv,using_cubic)
      
   
   return trcv

# EDA2 : 
# see /home/msok/Desktop/EDA2/logbook/20200128_EDA2_receiver_temperature.odt
def trcv_eda2( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too 
   print "INFO : using trcv_eda2 (same as EDA-1 + receiver temperature of the FEM in the SmartBox"

   x=[ 40, 50, 60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350 ]
   y=[ 13517.9, 13517.9,11947.6,11202.1,10830.2,10539,10205.4,10015.1,9683.08,9575.04,9396.75,9193.21,9046.74,8947.86,8820.11,8807.01,8596.01,8451,8367.9,8316.72,8230.98,8134.13,8023.28,7987.7,7761.26,7889.29,7796.44,7724.62,7741.32,7572.9,7504.12,7483.41 ]

   gain_eda2_lna = 100.00 # about 20dB 
   trcv_fem = 10000.00 / 100.00 # Approximately 10000 K / 20 dB of MWA LNA
   if freq_mhz >= 40 and freq_mhz <= 350 :
       trcv_fem_interpol = interp1d(x, y, kind='cubic')
       trcv_fem = trcv_fem_interpol( freq_mhz ) / gain_eda2_lna

   # at the moment return the same as for EDA-1 (or EDA-2)
   trcv_eda_lna = trcv_from_skymodel_with_err( freq_mhz , use_cubic )
   trcv_out = trcv_eda_lna + trcv_fem
   print "INFO : using trcv_eda2( %.2f MHz) = %.2f + %.2f  = %.2f [K]" % (freq_mhz,trcv_eda_lna,trcv_fem,trcv_out)
   
   return trcv_out


# EDA1 (as for the EDA1 paper) :
def trcv_eda1( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too
   trcv_out = trcv_from_skymodel_with_err( freq_mhz , use_cubic )
   
   print "INFO : using trcv_eda1( %.2f MHz) = %.2f [K]" % (freq_mhz,trcv_out)
   
   return trcv_out


def trcv_from_skymodel_with_err_cubic( freq_mhz ) :
    x=[55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0]
    y=[3722.55055596,1421.83480972,655.8033454,409.69506365,356.59552554,293.4326318,200.26593869,163.91940732,157.46121101,112.19530177,51.14849704,68.56693587,28.33318976,43.72596282,34.41633474,62.76422256]
    tlna_cubic = interp1d(x, y, kind='cubic')
    trcv = tlna_cubic(freq_mhz)
    
    return trcv


def trcv_lightcurve201612( freq_mhz ):
    T_rcv =  +53480.71172570309136062860488891601562500000000000000000 - 1455.76561026112926811038050800561904907226562500000000*freq_mhz \
             +12.16264789947674174186431628186255693435668945312500*freq_mhz*freq_mhz -0.00271098517125182300949171043669139180565252900124*freq_mhz*freq_mhz*freq_mhz \
             -0.00030428122418899071835451941581140999915078282356*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000046193435183846278872557794374642536894270961*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000959650510252193562322564293230395260358989*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000001953073405836920809493335023823663365930*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000028926757335617184960944065203237193611*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000055485231182041931661150443032917814*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz
             
    return T_rcv             
   
def trcv_fit_lightcurve_20160703_polfit( freq_mhz ):
   t_lna =  +36853.38919790466752601787447929382324218750000000000000  -1377.64458069924739902489818632602691650390625000000000 * freq_mhz +20.37196908756931890138730523176491260528564453125000 * freq_mhz*freq_mhz \
            -0.13515251806813355361924777753301896154880523681641 * freq_mhz*freq_mhz*freq_mhz +0.00022215700342832035230414278181854115246096625924 * freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            +0.00000169842384409759492976880849990362065682347747 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000464297999395697167928294570157440879221156 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            -0.00000000003336389317067374617426693716906293402102 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000000017598078793649799250229898970650520328 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            -0.00000000000000023477806740472191603243644100388291 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   return t_lna


def trcv_multi( freq_mhz , type, use_cubic=False ):
   t_rcv = 0.00
   if type.lower() == "eda2" or type.lower() == "trcv_eda2" :
      t_rcv = trcv_eda2( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "eda1" or type.lower() == "trcv_eda1" :
      t_rcv = trcv_eda1( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "aavs2" or type.lower() == "trcv_aavs2" :
      t_rcv = trcv_aavs2( freq_mhz, use_cubic=use_cubic )       
   else :
      print "ERROR : unknown receiver temperature type = %s -> CRITICAL ERROR -> cannot continue" % (type)
      exit(-1)
      # t_rcv = trcv_from_skymodel_with_err_cubic( freq_mhz )

   return t_rcv


if __name__ == "__main__":
    freq_mhz = 160.00
    if len(sys.argv) >= 1:
       freq_mhz=float( sys.argv[1] )

    type="eda2"
    if len(sys.argv) >= 2:
       type=sys.argv[2]

    t_rcv = trcv_multi( freq_mhz , type, False )
    print "%s : T_rcv (%.2f MHz) = %.2f [K]" % (type,freq_mhz,t_rcv)
       
       
    
    