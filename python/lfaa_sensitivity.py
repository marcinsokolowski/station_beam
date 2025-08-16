import lfaa_requirements
import sys
import math

# options :
from optparse import OptionParser,OptionGroup


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tSensitivity of the SKA-Low station\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('--n_polarisations','--n_pols','--npols',dest="n_polarisations",default=2, help="Number of polarisations [default %default]",metavar="int")
   parser.add_option('--lofar_rate',dest="lofar_rate",default=False,action="store_true", help="Use LOFAR rate [default %default]")
   parser.add_option('--lofar_rate_value',dest="lofar_rate_value",default=3, help="Use LOFAR rate value [default %default]",type="int")
   parser.add_option('--scaling_index',dest="scaling_index",default=-2.1, help="Use LOFAR rate value [default %default]",type="float")
   parser.add_option('--euclidean',dest="euclidean",default=False,action="store_true", help="Euclidean source count scaling [default %default]")
   parser.add_option('--verb',dest="verb",default=True,action="store_true", help="Verbosity level [default %default]")
   
   (options, args) = parser.parse_args(sys.argv[idx:])
   
   return (options, args)


# see equation 6 in Sensitivity DB paper (Sokolowski et al. (2021))
# sens_jy = sefd_station / math.sqrt( bw_hz * inttime_sec * options.n_polarisations )
def imaging_sensitivity( sefd_station, bandwidth_hz=30720000, inttime_sec=120, antnum=512, efficiency=1.00 ) :
    image_noise = -1000
    if antnum >= 2 :
       image_noise = (efficiency * sefd_station) / math.sqrt(bandwidth_hz * inttime_sec * antnum * (antnum - 1) * 0.5 ) # missing factor of 1/2 added from number of baselines, see Tools of Radio Astronomy page 228
    else :
       image_noise = (efficiency * sefd_station) / math.sqrt(bandwidth_hz * inttime_sec )
       
    
    return image_noise


def frb_rate( fluence_min, min_elev_deg=20, scaling_index=-2.1, verb=False ):
   # Shannon et al. 2018 
   # 37 +/- 8 /day/sky 
   # above threshold = 26 Jy ms (w/1.26 ms)^{-1/2}
   
   # min_elevation = 20 degrees is what is sensible for the EDA2/AAVS2 all-sky images:
   
   
   per_sky_per_day = 37*(pow( (fluence_min/26.00) , scaling_index ))
   per_sky_per_year = per_sky_per_day*365
   
   if verb :
      print("per_sky_per_day = %.4f FRBs / sky / day" % (per_sky_per_day))
   
   # multiply by the fraction of the sky corresponding to elevation range (min_elevation,90) [deg]
   min_elev_rad = min_elev_deg * (math.pi / 180.00)    
   sky_fraction = 0.5*(1.00 - math.sin( min_elev_rad ) )
#   print("DEBUG : min_elevation = %.2f [deg] -> fraction = %.5f" % (min_elev_deg,sky_fraction))
   per_sky_per_day = per_sky_per_day * sky_fraction
   per_sky_per_year = per_sky_per_year * sky_fraction
   
   return (per_sky_per_day,per_sky_per_year)

# from paper : https://ui.adsabs.harvard.edu/abs/2020arXiv201208348P/abstract
# combining our results with previous upper-limits on the all-sky FRB rate at 150 MHz, we find that there are 3-450 FRBs/sky/day above 50 Jy ms at 90% confidence.
# index=-2.1 as in Ryan Shannon et al. (2018) , Euclidean space is -1.5 = -3/2 
def lofar_frb_rate( fluence_min, min_elev_deg=20, lofar_rate=3, cut_off_fluence=50, verb=False, scaling_index=-2.1 ):

   per_sky_per_day =  lofar_rate*(pow( (float(fluence_min)/float(cut_off_fluence)) , scaling_index ))
   per_sky_per_year = per_sky_per_day*365
   
   if verb :
      print("per_sky_per_day = %.4f FRBs / sky / day" % (per_sky_per_day))
      

   # multiply by the fraction of the sky corresponding to elevation range (min_elevation,90) [deg]
   min_elev_rad = min_elev_deg * (math.pi / 180.00)
   sky_fraction = 0.5*(1.00 - math.sin( min_elev_rad ) )
   per_sky_per_day = per_sky_per_day * sky_fraction
   per_sky_per_year = per_sky_per_year * sky_fraction

   return (per_sky_per_day,per_sky_per_year)


if __name__ == "__main__":
    freq_mhz = 160.00
    if len(sys.argv) >= 1 :
       freq_mhz = float( sys.argv[1] )

    inttime=1000.00 # ms 
#    if len(sys.argv) >= 3 :
#       inttime = float( sys.argv[2] )
       
    bw_chan=(400.00/512.00) # *(32.00/27.00) # single channel :    
    n_chan=1
    if len(sys.argv) >= 3 :
       n_chan = int( sys.argv[2] )
    bw = bw_chan*n_chan   

    (options, args) = parse_options()
    
    if options.euclidean :
       options.scaling_index = -1.5


    print("###############################")
    print("PARAMETERS:")
    print("###############################")
    print("Frequency        = %.2f MHz" % (freq_mhz))
    print("Integration time = %.2f ms" % (inttime))
    print("N_channels       = %d -> %d x %.2f MHz = %.2f MHz" % (n_chan,n_chan,bw_chan,bw))
    print("N polarisations  = %d" % (options.n_polarisations))
    print("Use LOFAR Pastor-Marazuela (2020) values = %s [value = %d / day / sky]" % (options.lofar_rate,options.lofar_rate_value))
    print("Source count scaling index = %.4f (Euclidean flag = %s)" % (options.scaling_index,options.euclidean))
    print("###############################")
           
    
    aot_station = lfaa_requirements.lfaa_per_station( freq_mhz , interpolation_kind='cubic')
    sefd_station = lfaa_requirements.aot2sefd( aot_station )
    
    print("Frequency = %.2f MHz" % (freq_mhz))
    print("Station A/T = %.2f m^2/K" % (aot_station))
    print("Station SEFD = %.2f m^2/K" % (sefd_station))
 
    bw_hz = bw*1000000.00
    print("DEBUG : bw_hz = %.6f [Hz]" % (bw_hz))
    
    for n_sigma in (3,5,8,10) :
       outname = "sens_nchan%d_freq%.2fMHz_%.1fsigma.txt" % (n_chan,freq_mhz,n_sigma)
       out_f = open( outname , "w" )
       
       # (inttime_ms,limit_ms_Nsigma,per_sky_per_day,per_sky_per_year,sens_mjy,n_sigma)
       line = "# INTTIME[ms] Fluence_limit[Jy ms] FRBs/day FRBs/year Sens[mJy] N_sigma\n"
       out_f.write( line )
    
#       for inttime_ms in (1,10,50,100,150,200,300,500,762,1000,2000,2286,3000,4000,5000,10000) : 
       for inttime_ms in (1,10) :
          inttime_sec = inttime_ms/1000.00
              
          sens_jy = sefd_station / math.sqrt( bw_hz * inttime_sec * options.n_polarisations )
          sens_mjy = sens_jy*1000.00
          print("DEBUG : sens_jy = %.8f [Jy]" % (sens_jy))
       
#       n_sigma=3
#       limit_ms_3sigma = sens_jy*n_sigma*inttime_sec
#       limit_ms_3sigma2 = sefd_station * math.sqrt( inttime_sec / bw_hz ) * n_sigma             
#       print("\t%.1f ms : %.2f mJy -> 3sigma limit = %.2f [Jy sec] (vs. %.2f)" % (inttime_ms,sens_mjy,limit_ms_3sigma,limit_ms_3sigma2))
    

          limit_ms_Nsigma = sens_jy*n_sigma*inttime_ms
          
          if options.lofar_rate :
             (per_sky_per_day,per_sky_per_year) = lofar_frb_rate( limit_ms_Nsigma, lofar_rate=options.lofar_rate_value, scaling_index=options.scaling_index , verb=options.verb )
          else :
             (per_sky_per_day,per_sky_per_year) = frb_rate( limit_ms_Nsigma, scaling_index=options.scaling_index, verb=options.verb )
          
          print("\t%.1f ms : %.2f mJy -> %dsigma limit = %.2f [Jy msec] -> %.5f / day / sky or %.5f / year / sky" % (inttime_ms,sens_mjy,n_sigma,limit_ms_Nsigma,per_sky_per_day,per_sky_per_year))
          
          line = "%1.f %.2f %.5f %.5f %.2f %.1f\n" % (inttime_ms,limit_ms_Nsigma,per_sky_per_day,per_sky_per_year,sens_mjy,n_sigma)
          out_f.write( line )
         
    
#       n_sigma=10
#       limit_ms_3sigma = sens_jy*n_sigma*inttime_ms
#       print("\t%.1f ms : %.2f mJy -> %dsigma limit = %.2f [Jy msec]" % (inttime_ms,sens_mjy,n_sigma,limit_ms_3sigma))


       print("")    
       out_f.close()
    