#!/usr/bin/python
"""
primarybeammap.py --freq=202.24 --beamformer=0,0,0,1,3,3,3,3,6,6,6,6,8,9,9,9 --datetimestring=20110926210616

main task is:
make_primarybeammap()

This is the script interface to the functions and modules defined in MWA_Tools/src/primarybeamap.py

"""
import pdb

# from mwapy import ephem_utils
# from mwapy.pb import primary_beam
import sys
import numpy,math
import os,time
from optparse import OptionParser
import ephem
import logging
import matplotlib
if not 'matplotlib.backends' in sys.modules:
    matplotlib.use('agg')
import matplotlib.pyplot as pylab
import astropy.io.fits as pyfits
#from mwapy.pb.primarybeammap import *
from primarybeammap_local import *
#import primarybeammap
import errno
import tpm_mapping

# receiver temperature module :
import station_trcv

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise  

def read_antenna_list( filename, start_column=0 ):
   file = open( filename , 'r' )
   data = file.readlines()
   
   x_arr = []
   y_arr = []
   z_arr = []
   count = 0
   for line in data : 
       words = line.split() # was ' ' , but when none is provided -> all white-space characters are ok !
       print "DEBUG : words = %s (len = %d)" % (words,len(words))
       
       if line[0] == '#' :
          continue

       if line[0] != "#" :
           x = float(words[0+start_column])
           y = float(words[1+start_column])
           z = 0            
           if len(words) >= 3 :
              z = float(words[2+start_column])
      
           x_arr.append(x)
           y_arr.append(y)
           z_arr.append(z)
           count = count + 1

   file.close()
   
   print "Read %d / %d / %d of x, y, z positions from file %s" % (len(x_arr),len(y_arr),len(z_arr),filename)
           
   return (x_arr,y_arr,z_arr,count)         


     

def main():

    usage="Usage: %prog [options]\n"
    usage+="\tCreates an image of the 408 MHz sky (annoted with sources) that includes contours for the MWA primary beam\n"
    usage+="\tThe beam is monochromatic, and is the sum of the XX and YY beams\n"
    usage+="\tThe date/time (UT) and beamformer delays must be specified\n"
    usage+="\tBeamformer delays should be separated by commas\n"
    usage+="\tFrequency is in MHz, or a coarse channel number (can also be comma-separated list)\n"
    usage+="\tDefault is to plot centered on RA=0, but if -r/--racenter, will center on LST\n"
    usage+="\tContours will be plotted at %s of the peak\n" % contourlevels
    usage+="\tExample:\tpython primarybeammap.py -c 98 --beamformer=1,0,0,0,3,3,3,3,6,6,6,6,9,9,9,8 --datetimestring=20110926211840\n\n"
    
    parser = OptionParser(usage=usage)
    parser.add_option('-d','--datetimestring',dest="datetimestring",default=None,
                      help="Compute for <DATETIMESTRING> (YYYYMMDDhhmmss)",
                      metavar="DATETIMESTRING")
    parser.add_option('-c','--channel',dest='channel',default=None,
                      help='Center channel(s) of observation')
    parser.add_option('-f','--frequency',dest='frequency',default="150",
                      help='Center frequency(s) of observation [MHz]')
    parser.add_option('-b','--beamformer',dest='delays',default="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
                      help='16 beamformer delays separated by commas')
    parser.add_option('-D','--date',dest='date',default=None,
                      help='UT Date')
    parser.add_option('-t','--time',dest='time',default=None,
                      help='UT Time')
    parser.add_option('-g','--gps',dest='gps',default=0,help='GPS time')
    parser.add_option('-m','--model',dest='model',default='analytic',
                      help='beam model: analytic, advanced, full_EE, full_EE_AAVS05')       
    parser.add_option('-p','--plottype',dest='plottype',default='beamsky',
                      help='Type of plot: all, beam, sky, beamsky, beamsky_scaled')
    parser.add_option('--title',dest='title',default=None,
                      help='Plot title')
    parser.add_option('-e','--ext',dest='extension',default='png',
                      help='Plot extension [default=%default]')
    parser.add_option('-r','--racenter',action="store_true",dest="center",default=False,
                      help="Center on LST?")
#    parser.add_option('-s','--sunline',dest="sunline",default="1",choices=['0','1'],
#                      help="Plot sun [default=%default]")
    parser.add_option('--tle',dest='tle',default=None,
                      help='Satellite TLE file')
    parser.add_option('--duration',dest='duration',default=300,type=int,
                      help='Duration for plotting satellite track')
    parser.add_option('--size',dest='size',default=1000,type=int,
                      help='Resolution of created beam file')
    parser.add_option('--dir',dest='dir',default=None,help='output directory')

    parser.add_option('-v','--verbose',action="store_true",dest="verbose",default=False,
                      help="Increase verbosity of output")

    parser.add_option('-u','--unixtime',dest='unixtime_start',default=1471305600,help='Unixtime start',type=int)
    parser.add_option('-i','--interval',dest='interval',default=86400,help='Unixtime start',type=int)
    parser.add_option('-s','--step',dest='step',default=1800,help='Step in seconds',type=int)
    parser.add_option('-o','--outfile',dest='out_filename',default='out.txt',help='Output filename')
    parser.add_option('-z','--za',dest='pointing_za_deg',default=0.00,help='Pointing za [deg]',type=float)
    parser.add_option('-a','--az',dest='pointing_az_deg',default=0.00,help='Pointing az [deg]',type=float)

    # polarisations :
    parser.add_option('--pols','--polarisations',dest='pols',default='XX',help='Polarisations [default %default]')

    # haslam map :
    parser.add_option('--haslam','--radio_image',dest='radio_image',default='radio408.RaDec.fits',help='Sky model map [default %default]')

    
    # sun:
    parser.add_option('--add_sun',action="store_true",dest="add_sun",default=False,help="Add sun")
    parser.add_option('--sun_ra',dest="sun_ra",default=290.46450056,help="Sun RA [deg] default = %default [deg]",type="float")
    parser.add_option('--sun_dec',dest="sun_dec",default=-26.70331900,help="Sun Dec [deg] default = %default [deg]",type="float")
#    parser.add_option('--sun_scale_size',dest="sun_scale_size",default=1.7,help="Scale size of the Sun to be more like radio image (not optical)",type="float")
    parser.add_option('--sun_scale_size',dest="sun_scale_size",default=1.0,help="Scale size of the Sun to be more like radio image (not optical)",type="float")
    
    # gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00
    parser.add_option('--gain_sigma_db',dest='gain_sigma_db',default=0,help="RMS of the gain of each dipole",type=float)
    parser.add_option('--gain_sigma_ph',dest='gain_sigma_ph',default=0,help="RMS of the phase at 160 MHz of each dipole",type=float)
            
    # receiver noise :
    parser.add_option('--use_trcv',action="store_true",dest="use_trcv",default=False,help="Use T_rcv [default %s]")
    parser.add_option('--T_rcv','--t_rcv',dest='T_rcv',default=0.00,help="Receiver noise [default %]",type=float)
    parser.add_option('--rcv_noise_coupling',dest='rcv_noise_coupling',default=0.00,help="Receiver noise coupling [default %]",type=float)
    parser.add_option('--trcv_type',dest='trcv_type',default=None,help="Source of T_rcv possible values = lightcurve201612 (default), lightcurve_20160703_cubic, lightcurve_20160703_polfit, budi, data_vs_model_201612, trcv_angelica_data_vs_time, trcv_angelica_data_vs_time_powerlawfit, trcv_from_skymodel_with_err, trcv_aavs2, trcv_aavs2_vs_za_deg, trcv_eda2, trcv_eda1")
 

    # only include a sub-array (list of tpm-s is more for the AAVS1 but also EDA-2 ):
    parser.add_option('--ant_list','--antenna_list',dest='antenna_list',default=None, help='List of antennas to be simulated (as a sub-array) [default %]')
    parser.add_option('--tpm_list',dest='tpm_list',default=None, help='List of TPMs to be simulated (as a sub-array [default %]')
    parser.add_option('--ant_list_file','--antenna_list_file',dest='antenna_list_file',default=None, help='File with a list of antenna positions to be simulated [default %]')
    parser.add_option('--start_column','--antlist_start_column',dest='antlist_start_column',default=0, help='[default %]',type="int")

    
    # feko file :
    parser.add_option('--feko_file','--feko',dest='feko_file',default=None, help='FEKO File [default %]')
    parser.add_option('--interpolate',action="store_true",dest="interpolate",default=False, help="Interpolate beam [default %]")
    parser.add_option('--interpolate_to_nearest',action="store_true",dest="interpolate_to_nearest",default=False, help="Interpolate beam to nearest value [default %]")

    # Beam models using fits file with a single, isolated, dipole pattern :
    # use beam fits files in ~/aavs-calibration/
    # use_beam_fits=options.use_beam_fits, station_name=options.station_name
    parser.add_option('--use_beam_fits',action="store_true",dest="use_beam_fits",default=False,help="Use fits beam files from ~/aavs-calbration/BeamModels/ [default %default]")
    parser.add_option('--beam_fits_station','--station_name',dest="station_name",default="EDA",help="Station name to find beam fits files in ~/aavs-calbration/BeamModels/ [default %default]")

    
    tpm_list=None
    xpos_list = None
    ypos_list = None
    zpos_list = None
    (options, args) = parser.parse_args()
    datetimestring=options.datetimestring
    
    if options.trcv_type is not None :
       if not options.use_trcv :
          options.use_trcv = True
          print "WARNING : forcing receiver temperature flag use_trcv to True because trcv_type = %s (!= None)" % (options.trcv_type)

    set_sun( options.add_sun, options.sun_ra, options.sun_dec, scale_size=options.sun_scale_size )

    if options.tpm_list is not None :
        tpm_list = options.tpm_list.split(",")
        tpm_list = map(int,tpm_list)
        
        eda_station = tpm_mapping.Station()
        (xpos_list,ypos_list) = eda_station.get_xy_pos( tpm_list[0] )
        print "Antenna positions :"
        print "\txs = %s" % (xpos_list)
        print "\tys = %s" % (ypos_list)
        

    
    print "##################################################"
    print "PARAMETERS:"
    print "##################################################"
    print "Interval    = %d [sec]" % (options.interval)
    print "Add sun     = %s at (ra,dec) = (%.4f,%.4f) [deg] , size scaled by x %.4f " % (options.add_sun,options.sun_ra,options.sun_dec,options.sun_scale_size)
    print "Beamoforming errors RMS(gain) = %.4f dB , RMS(phase) = %.4f degree" % (options.gain_sigma_db,options.gain_sigma_ph)
    print "Use T_rcv   = %s" % (options.use_trcv)
    print "T_rcv type  = %s" % (options.trcv_type)
    print "Receiver noise = %.2f [K]" % (options.T_rcv)
    print "Receiver noise coupling = %.2f" % (options.rcv_noise_coupling)
    print "TPM list     = %s -> %s" % (options.tpm_list,tpm_list)
    if options.tpm_list is not None :
        print"\txs = %s" % (xpos_list)
        print"\tys = %s" % (ypos_list)
    print "Antenna list = %s" % (options.antenna_list)
    print "Antenna list file = %s" % (options.antenna_list_file)
    print "Antenna list start column = %d" % (options.antlist_start_column)
    print "FEKO file         = %s" % (options.feko_file)
    print "Sky image resolution = %d x %d" % (options.size,options.size)
    print "Interpolate       = %s (to nearest = %s)" % (options.interpolate,options.interpolate_to_nearest)
    print "Polarisations     = %s" % (options.pols)
    print "##################################################"
            
    
    if options.dir is not None:
       mkdir_p(options.dir)

    
    if options.frequency is not None:
        if (',' in options.frequency):
            try:
                frequency=map(float,options.frequency.split(','))
            except ValueError:
                logger.error("Could not parse frequency %s\n" % options.frequency)
                sys.exit(1)
        else:
            try:
                frequency=float(options.frequency)
            except ValueError:
                logger.error("Could not parse frequency %s\n" % options.frequency)
                sys.exit(1)
    else:
        frequency=options.frequency
    if options.channel is not None:
        if (',' in options.channel):
            try:
                channel=map(float,options.channel.split(','))
            except ValueError:
                logger.error("Could not parse channel %s\n" % options.channel)
                sys.exit(1)
        else:
            try:
                channel=float(options.channel)
            except ValueError:
                logger.error("Could not parse channel %s\n" % options.channel)
                sys.exit(1)
    else:
        channel=options.channel
    if options.delays is not None:
        try:
            if (',' in options.delays):
                delays=map(int,options.delays.split(','))
            else:
                delays=16*[int(options.delays)]
        except:
            logger.error("Could not parse beamformer delays %s\n" % options.delays)
            sys.exit(1)
    else:
        delays=options.delays
    extension=options.extension
    plottype=options.plottype
    model=options.model
    if model not in ['analytic','advanced','full_EE', 'full_EE_AAVS05']:
        logger.error("Model %s not found\n" % model)
        sys.exit(1)   
    if plottype.lower() not in ['all', 'beam','sky','beamsky','beamsky_scaled','none']:
        logger.error("Plot type %s not found\n" % plottype)
        sys.exit(1)                   
    gpsstring=options.gps
    gps=int(gpsstring)
               
    if (len(delays)<16):
        logger.error("Must supply 1 or 16 delays\n")
        sys.exit(1)
    if (frequency is None):
        if (channel is not None):
            if (isinstance(channel,list)):
                frequency=list(1.28*numpy.array(channel)) # multiplication by 1e6 is done later at line Convert to Hz
            else:
                frequency=1.28*channel # multiplication by 1e6 is done later at line Convert to Hz
    if frequency is None:
        logger.error("Must supply frequency or channel\n")
        sys.exit(1)
    if (isinstance(frequency,int) or isinstance(frequency,float)):
        frequency=[frequency]
    frequency=np.array(frequency)*1e6 #Convert to Hz        

    out_file=open(options.out_filename,"w")
    out_file.close()

    if options.antenna_list_file is not None :
        print "Reading antenna positions from file : %s" % (options.antenna_list_file)
        
        (xpos_list,ypos_list,zpos_list,ant_count) = read_antenna_list( options.antenna_list_file, start_column=options.antlist_start_column )
        print "Read %d antennas from file %s, list of antennas:" % (ant_count,options.antenna_list_file)        
        print "\txs = %s" % (xpos_list)
        print
        print "\tys = %s" % (ypos_list)
        print
        print "\tzs = %s" % (zpos_list)
             
    # calculate receiver temperature :
    T_rcv = options.T_rcv
    if options.trcv_type is not None :
       T_rcv = station_trcv.trcv_multi( frequency/1e6 , options.trcv_type )
       print "Receiver temperature T_rcv( %.2f MHz, type=%s ) = %.4f [Kelvin]" % ((frequency/1e6),type,T_rcv)
       

    beams = None
    for uxtime in range(options.unixtime_start,options.unixtime_start+options.interval,options.step) :
        gps = uxtime - 315964800
        print "unixtime = %d -> gpstime = %d" % (uxtime,gps)
        
        for freq in frequency:
           print 'frequency', freq
           (beamsky_sum_XX,beam_sum_XX,Tant_XX,beam_dOMEGA_sum_XX,beamsky_sum_YY,beam_sum_YY,Tant_YY,beam_dOMEGA_sum_YY,beams) = make_primarybeammap( gps, delays, freq, beams=beams, model=model, plottype=plottype, extension=extension, resolution=options.size, directory=options.dir, out_filename=options.out_filename,
                                       pointing_za_deg=options.pointing_za_deg, pointing_az_deg=options.pointing_az_deg, gain_sigma_ph_160mhz=options.gain_sigma_ph, gain_sigma_dB=options.gain_sigma_db,
                                       radio_image=options.radio_image,
                                       T_rcv=T_rcv, rcv_noise_coupling=options.rcv_noise_coupling, use_trcv=options.use_trcv, xpos=xpos_list, ypos=ypos_list, zpos=zpos_list,
                                       feko_file=options.feko_file, interpolate=options.interpolate, interpolate_to_nearest=options.interpolate_to_nearest, pol_list_string=options.pols,
                                       use_beam_fits=options.use_beam_fits, station_name=options.station_name )
                                       
           # to keep the normal EDA version as it was - beam regenerated in every call
           # TODO : once tested it can be changed 
#           if options.feko_file is None :
#               beams = None                                       
        
        
                                       
#           if (result is not None):
#               print "Wrote %s" % result
    
#    alpha=-2.55
#    rX,rY=get_skytemp(datetimestring, delays, frequency,alpha=alpha)
#    print 'Tsky_X for alpha %s: %s'%(alpha,rX)
#    print 'Tsky_Y for alpha %s: %s'%(alpha,rY)

if __name__ == "__main__":
    main()
