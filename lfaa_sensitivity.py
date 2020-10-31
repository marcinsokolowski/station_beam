import lfaa_requirements
import sys
import math

if __name__ == "__main__":
    freq_mhz = 160.00
    if len(sys.argv) >= 1 :
       freq_mhz = float( sys.argv[1] )

    inttime=1000.00 # ms 
    if len(sys.argv) >= 3 :
       inttime = float( sys.argv[2] )
       
    bw_chan=(400.00/512.00)*(32.00/27.00) # single channel :    
    n_chan=1
    if len(sys.argv) >= 4 :
       n_chan = int( sys.argv[3] )
    bw = bw_chan*n_chan   

    
    print("###############################")
    print("PARAMETERS:")
    print("###############################")
    print("Frequency = %.2f MHz" % (freq_mhz))
    print("Integration time = %.2f ms" % (inttime))
    print("N_channels = %d -> %d x %.2f MHz = %.2f MHz" % (n_chan,n_chan,bw_chan,bw))
    print("###############################")

           
    
    aot_station = lfaa_requirements.lfaa_per_station( freq_mhz , interpolation_kind='cubic')
    sefd_station = lfaa_requirements.aot2sefd( aot_station )
    
    print("Frequency = %.2f MHz" % (freq_mhz))
    print("Station A/T = %.2f m^2/K" % (aot_station))
    print("Station SEFD = %.2f m^2/K" % (sefd_station))
 
    bw_hz = bw*1000000.00
    for inttime_ms in (10,100,200,500,1000,2000) : 
       inttime_sec = inttime_ms/1000.00
              
       sens_jy = sefd_station / math.sqrt( bw_hz * inttime_sec )
       sens_mjy = sens_jy*1000.00
       
       n_sigma=3
#       limit_ms_3sigma = sens_jy*n_sigma*inttime_sec
#       limit_ms_3sigma2 = sefd_station * math.sqrt( inttime_sec / bw_hz ) * n_sigma             
#       print("\t%.1f ms : %.2f mJy -> 3sigma limit = %.2f [Jy sec] (vs. %.2f)" % (inttime_ms,sens_mjy,limit_ms_3sigma,limit_ms_3sigma2))
    
       n_sigma=3
       limit_ms_3sigma = sens_jy*n_sigma*inttime_ms
       print("\t%.1f ms : %.2f mJy -> %dsigma limit = %.2f [Jy msec]" % (inttime_ms,sens_mjy,n_sigma,limit_ms_3sigma))
    
       n_sigma=10
       limit_ms_3sigma = sens_jy*n_sigma*inttime_ms
       print("\t%.1f ms : %.2f mJy -> %dsigma limit = %.2f [Jy msec]" % (inttime_ms,sens_mjy,n_sigma,limit_ms_3sigma))
       print("")
    
