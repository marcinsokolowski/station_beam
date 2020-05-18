from __future__ import print_function
# from pylab import *
import pdb

import numpy
import os,sys
import math
from string import Template
import beam_tools

# import pyfits
import astropy.io.fits as pyfits

current_fits_filename = None
current_fits_beam     = None 
azim_map              = None
za_map                = None


def save_fits( data , out_fits_name ) :
   hdu = pyfits.PrimaryHDU()
   hdu.data = data
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , clobber=True )

   print("Saved fits file %s" % (out_fits_name))   

# dist_za = numpy.sin( za_deg*(math.pi/180.00) )*(npix/2)
# x_pixel = ( npix/2 - dist_za*numpy.sin( azim_deg*(math.pi/180.00) ) ) 
# y_pixel = ( npix/2 + dist_za*numpy.cos( azim_deg*(math.pi/180.00) ) )   
# print "x_pixel : count = %d" % (len(x_pixel))

# def find_value_fast( az, za, map_az, map_za, x_map, y_map, map_beam, verb=False ) :

def find_value_fast( az, za, x_start, y_start, map_az, map_za, map_beam, radius=20, verb=False ) :
   x_size = map_az.shape[0]
   y_size = map_az.shape[1]

   min_distance = 1000000.000
   beam_value   = -1e6


   ra1  = az
   dec1 = (math.pi/2 - za)

   sin_dec1 = math.sin(dec1)
   cos_dec1 = math.cos(dec1)
   
   print("\tfind_value_fast : around (%d,%d) in range (%d,%d) - (%d,%d)" % ( x_start,y_start, max(x_start-radius,0),max(y_start-radius,0),min(x_start+radius,x_size),min(y_start+radius,y_size)))
   best_x = -1
   best_y = -1

   
   for y in range( max(y_start-radius,0) , min(y_start+radius,y_size) ) :
      for x in range( max(x_start-radius,0) , min(x_start+radius,x_size) ) :
         if not numpy.isnan( map_az[y,x] ) and not numpy.isnan( map_za[y,x] ) :
             ra2  = map_az[y,x]
             dec2 = (math.pi/2 - map_za[y,x])

             cos_value = sin_dec1*math.sin(dec2) + cos_dec1*math.cos(dec2)*math.cos( ra1 - ra2 )
             ang_distance = math.acos( cos_value )

             if verb :
                print("\t(%d,%d) -> ang_distance = %.4f [arcsec]" % (x,y,ang_distance*(180.00/math.pi)*3600))

             if ang_distance < min_distance :
                min_distance = ang_distance 
                beam_value = map_beam[y,x]
                best_x = x
                best_y = y


                if verb :
                   print("New minimum at (%d,%d) - distance = %.2f [arcsec] , beam_value = %.2f" % (x,y,min_distance*(180.00/math.pi)*3600,beam_value))

   print("\tfind_value_fast : best (x,y) = (%d,%d) , min_distance = %.2f [arcsec]" % (best_x,best_y,min_distance*(180.00/math.pi)*3600))

   return (beam_value,min_distance*(180.00/math.pi)*3600)
   

def find_value( az, za, map_az, map_za, map_beam, verb=False ) :
   x_size = map_az.shape[0]
   y_size = map_az.shape[1]
   
   min_distance = 1000000.000
   beam_value   = -1e6

   ra1  = az
   dec1 = (math.pi/2 - za)
   
   sin_dec1 = math.sin(dec1)
   cos_dec1 = math.cos(dec1)
   
   best_x = -1
   best_y = -1
   
   for y in range(0,y_size) :
      for x in range(0,x_size ) :
         if not numpy.isnan( map_az[y,x] ) and not numpy.isnan( map_za[y,x] ) :
            ra2  = map_az[y,x]
            dec2 = (math.pi/2 - map_za[y,x])
         
            cos_value = sin_dec1*math.sin(dec2) + cos_dec1*math.cos(dec2)*math.cos( ra1 - ra2 )
            ang_distance = math.acos( cos_value )
         
            if ang_distance < min_distance :
               min_distance = ang_distance 
               beam_value = map_beam[y,x]
               best_x = x
               best_y = y
            
               if verb :
                  print("New minimum at (%d,%d) - distance = %.2f [arcsec] , beam_value = %.2f" % (x,y,min_distance*(180.00/math.pi)*3600,beam_value))

   print("\tfind_value : best (x,y) = (%d,%d) , min_distance = %.2f [arcsec]" % (best_x,best_y,min_distance*(180.00/math.pi)*3600))

   return (beam_value,min_distance*(180.00/math.pi)*3600)
            
         
def remap_beam( fits_file , step=1, do_test=False, radius=20 ) :
   beam   = pyfits.open( fits_file )
   x_size = beam[0].header['NAXIS1']  
   y_size = beam[0].header['NAXIS2']
   npix = x_size 
   print("Read file %s with size (%d,%d)" % (fits_file,x_size,y_size))
   
   (az_sin,za_sin,x_sin,y_sin,z_sin,d_sin) = beam_tools.makeAZZA( x_size, 'SIN' ,azim_from_north=True, return_all=True )
   (az_zea,za_zea,x_zea,y_zea,z_zea,d_zea) = beam_tools.makeAZZA( x_size, 'ZEA' ,azim_from_north=True, return_all=True )
   
   print("Test orientation azim SIN :")
   print("         %.2f            " % (az_sin[npix-1,npix/2]))                  # WARNING data[y,x] !!!
   print("%.2f                 %.2f" % (az_sin[npix/2,0],az_sin[npix/2,npix-1])) # WARNING data[y,x] !!!   
   print("         %.2f            " % (az_sin[0,npix/2]))                       # WARNING data[y,x] !!! 

   print("Test orientation xy SIN :")
   print("         (%d,%d)            "   % (x_zea[npix/2,npix-1],y_zea[npix/2,npix-1]))
   print("(%d,%d)                (%d,%d)" % (x_zea[0,npix/2],y_zea[0,npix/2],x_zea[npix-1,npix/2],y_zea[npix-1,npix/2]))
   print("         (%d,%d)            "   % (x_zea[npix/2,0],y_zea[npix/2,0]))
   
#   dist_za_sin = numpy.sin( za_sin )*(npix/2)
#   x_pixel = ( npix/2 - dist_za_sin*numpy.sin( az_sin ) ) 
#   y_pixel = ( npix/2 + dist_za_sin*numpy.cos( az_sin ) )   

   
   save_fits(az_sin,"test_az_sin.fits")
   save_fits(za_sin,"test_za_sin.fits")
   
   save_fits(az_zea,"test_az_zea.fits")
   save_fits(za_zea,"test_za_zea.fits")
   
   out_file = fits_file.replace(".fits","_test.fits" )
   save_fits( beam[0].data, out_file )
      
   out_beam = numpy.zeros( beam[0].data.shape )
   for y in range(0,y_size,step):
      # print "i = %d" % (i)
      for x in range(0,x_size): 
         print("TEST (x,y) = (%d,%d)" % (x,y))
         az_zea_ij = az_zea[y,x] # this order of x,y is ok (see makeAZZA)
         za_zea_ij = za_zea[y,x] # this order of x,y is ok (see makeAZZA)
#         radius_zea = math.sqrt(1.00-math.cos(za_zea_ij))
         
         if not numpy.isnan( za_zea_ij ) :
             d_sin     = math.sin( za_zea_ij )
             y_zea_ij  = y_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
             x_zea_ij  = x_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
         
#             x_sin_ij = x_zea_ij * (d_sin/d_zea[i,j])
#             y_sin_ij = y_zea_ij * (d_sin/d_zea[i,j])
             d_sin = (npix/2)*math.sin( za_zea_ij )
#             r = ( radius_sin / radius_zea )
             y_sin_ij = d_sin*math.cos( az_zea_ij ) + (npix/2)
             x_sin_ij = d_sin*math.sin( az_zea_ij ) + (npix/2)
             
#             print "\t(%d,%d) : ZEA : (az,za) = (%.2f,%.2f) [deg], (%d,%d) r_zea = %.4f vs. SIN : (%.1f,%.1f) r_sin = %.4f -> (r_sin/r_zea) = %.8f" % (i,j,(az_zea_ij*(180.00/math.pi)),(za_zea_ij*(180.00/math.pi)),x_zea_ij,y_zea_ij,radius_zea,x_sin_ij,y_sin_ij,radius_sin,r)

# ??? WHY ZEA (x_zea,y_zea) = (%d,%d) is still other way around ???
         
             if do_test :     
                 print("\t(%d,%d) : ZEA (x_zea,y_zea) = (%d,%d) , (az,za) = (%.2f,%.2f) [deg] -> SIN : (x_sin,y_sin) = (%d,%d)" % (x,y,x_zea_ij,y_zea_ij,az_zea_ij,za_zea_ij,int(x_sin_ij),int(y_sin_ij)))

             # find given (az,za) from ZEA map in SIN map and return beam value (in SIN mapping) :
             beam_value = numpy.NaN
             if do_test :
                (beam_value, min_distance_arcsec ) = find_value( az_zea_ij, za_zea_ij, az_sin, za_sin, beam[0].data )                  
                (beam_value_fast, min_distance_arcsec_fast ) = find_value_fast( az_zea_ij, za_zea_ij, int(x_sin_ij), int(y_sin_ij), az_sin, za_sin, beam[0].data, radius=radius )
             
                if math.fabs(beam_value-beam_value_fast) > 0.00001 :
                   print("\n\n\t!!!! ERROR in code at (x,y) = (%d,%d) %.6f != %.6f !!!!\n\n" % (x,y,beam_value,beam_value_fast))
             else :
                # final function :
                (beam_value, min_distance_arcsec_fast ) = find_value_fast( az_zea_ij, za_zea_ij, int(x_sin_ij), int(y_sin_ij), az_sin, za_sin, beam[0].data, radius=radius )
                
             
             out_beam[y,x] = beam_value
         
         

   out_file = fits_file.replace(".fits","_zea.fits" )    
   save_fits( out_beam, out_file )            
         
         


def read_beam_fits( frequency_mhz, polarisation="X", station_name="EDA", simulation_path="$HOME/aavs-calibration/BeamModels/" , postfix="zea" ) :
   global current_fits_filename
   global current_fits_beam
   global azim_map
   global za_map

   polarisation_string = polarisation[0]
   postfix_full = ""
   if postfix is not None :
      postfix_full = "_" + postfix 
   beam_file_name_template = "%s/%s/%s/%s_%spol_ortho_%03d%s.fits" % (simulation_path,station_name,postfix.upper(),station_name,polarisation_string,int(frequency_mhz),postfix_full)
   beam_file_name = Template( beam_file_name_template ).substitute(os.environ)
   
      
   if current_fits_filename is None or current_fits_beam is None or beam_file_name != current_fits_filename :
      print("Reading fits file %s ..." % (beam_file_name))
           
      current_fits_beam = pyfits.open( beam_file_name )
      x_size = current_fits_beam[0].header['NAXIS1']  
      y_size = current_fits_beam[0].header['NAXIS2']

      print() 
      print('# \tRead fits file %s of size %d x %d' % (beam_file_name,x_size,y_size))
      current_fits_filename = beam_file_name
   else :
      print("Fits file name = %s is the same as previous %s -> not updating cache" % (beam_file_name,current_fits_filename))
      
      
#   x_size = current_fits_beam[0].header['NAXIS1']  
#   y_size = current_fits_beam[0].header['NAXIS2']
#   if azim_map is None or za_map is None :
#      (azim_map,za_map) = makeAZZA( x_size )
#      print "Generated local AZZA map for size %d x %d pixels" % (x_size,x_size)
  

   return (current_fits_beam,current_fits_filename)

# individual_antenna_beams_path="~/aavs-calibration/
# EDA_Xpol_ortho_169.fits
def get_fits_beam( azim_deg, za_deg, frequency_mhz, polarisation="X", station_name="EDA", simulation_path="$HOME/aavs-calibration/BeamModels/") :
   print("requestion beam model for azza map of size (%d x %d) pixels" % (azim_deg.shape[0],azim_deg.shape[1]))

   (current_fits_beam,beam_file_name) = read_beam_fits( frequency_mhz, polarisation, station_name, simulation_path )
   

   x_size = int( current_fits_beam[0].header['NAXIS1'] )
   y_size = int( current_fits_beam[0].header['NAXIS2'] )
   npix = x_size
   x_center = x_size / 2
   y_center = y_size / 2 
   print('INFO : Current fits file %s size = %d x %d -> center = %d x %d' % (beam_file_name,x_size,y_size,x_center,y_center))   

   # assuming that ds9 shows :   
   #       NORTH
   # EAST        WEST
   #       SOUTH 
   # FITS files are in orthographic projection :
   # d = np.sqrt(x*x + y*y)/(npix/2)
#   z = numpy.linspace(-npix/2.0,npix/2.0,num=npix) # was dtype=np.float32  
#   x = numpy.empty((npix,npix),dtype=numpy.float32)
#   y = numpy.empty((npix,npix),dtype=numpy.float32)
#   for i in range(npix):
#        y[i,0:] = z
#        x[0:,i] = z
#   d = numpy.sqrt(x*x + y*y)/(npix/2)

   # only select pixels above horizon
#   t = (d <= 1.0)
#   za = numpy.zeros((npix,npix),dtype=numpy.float32)*numpy.NaN
#   za[t] = numpy.arcsin(d[t])
#   print 'Using slant orthographic projection'
#   az = numpy.arctan2(y,x)
#   az = az+math.pi #0 to 2pi
#   az = 2*math.pi-az #Change to clockwise from top (when origin is in top-left)


   # see /home/msok/ska/aavs/aavs0.5/trunk/simulations/FEKO/beam_models/MWA_EE/MWAtools_pb/eda/beam_tools.py
   # get x,y from azim,za ... - other way around
   dist_za = numpy.sin( za_deg*(math.pi/180.00) )*(npix/2)
   x_pixel = ( npix/2 + dist_za*numpy.cos( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
   y_pixel = ( npix/2 + dist_za*numpy.sin( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
   print("x_pixel : count = %d" % (len(x_pixel)))

   # single pixel :
   beam_value = current_fits_beam[0].data[int(y_pixel),int(x_pixel)]
   print("Beam( %.4f deg, %.4f deg ) = %.8f at pixel (%d,%d)" % (azim_deg,za_deg,beam_value,x_pixel,y_pixel))
   return (beam_value)


   # any better way to check if it is an array ???
#   is_array = True
#   try :
#      print "len = %d" % (len(x_pixel))
#   except :
#      is_array = False
#   
#   if is_array : 
#       # multiple pixels :
#       x_pixel_lin = x_pixel.reshape( npix*npix )
#       y_pixel_lin = y_pixel.reshape( npix*npix )
#       x_pixel_lin_nonan = []
#       y_pixel_lin_nonan = []
#
#       for i in range(0,len(x_pixel_lin)) :
#           if not numpy.isnan( x_pixel_lin[i] ) and not numpy.isnan( y_pixel_lin[i] ) :
#               x_pixel_lin_nonan.append( x_pixel_lin[i] )
#               y_pixel_lin_nonan.append( y_pixel_lin[i] )
#   
#       x_pixel_lin_nonan = numpy.array( x_pixel_lin_nonan )
#       y_pixel_lin_nonan = numpy.array( y_pixel_lin_nonan )
#       cnt_nonan = len( x_pixel_lin_nonan )
#      
#       beam_values = numpy.zeros( (npix,npix) )
###   beam_values = current_fits_beam[0].data[y_pixel.astype(int),x_pixel.astype(int)]
##   for i in range(0,cnt_nonan) :
##      y_int = int(y_pixel_lin_nonan[i])
##      for j in range(0,cnt_nonan) :
##          x_int = int(x_pixel_lin_nonan[i])
###          print "(%d,%d) -> (%d,%d) = %.4f" % (i,j,x_int,y_int,current_fits_beam[0].data[y_int,x_int])
##          beam_values[x_int,y_int] = current_fits_beam[0].data[x_int,y_int] # or x,y ???
#       beam_values[x_pixel_lin_nonan.astype(int),y_pixel_lin_nonan.astype(int)] = current_fits_beam[0].data[x_pixel_lin_nonan.astype(int),y_pixel_lin_nonan.astype(int)] 
#       return (beam_values,current_fits_beam,x_pixel.astype(int),y_pixel.astype(int))
#   else :
#       # single pixel :
#       beam_value = current_fits_beam[0].data[int(y_pixel),int(x_pixel)]
#       print "Beam( %.4f deg, %.4f deg ) = %.8f at pixel (%d,%d)" % (azim_deg,za_deg,beam_value,x_pixel,y_pixel)
#       return (beam_value)

global_fits_counter  = 0

def makeAZZA(npix=256,projection='ZEA'):
    """
    Make azimuth and zenith angle arrays for a square image of side npix
    Projection is sine, all-sky
    Returns (az,za). Angles are in radian.
    """
    # build az and za arrays
    z = numpy.arange(npix,dtype=numpy.float32)-npix/2
    x = numpy.empty((npix,npix),dtype=numpy.float32)
    y = numpy.empty((npix,npix),dtype=numpy.float32)
    dOMEGA = numpy.empty((npix,npix),dtype=numpy.float32)
    
    for i in range(npix):
        y[i,0:] = z
        x[0:,i] = z
    d = numpy.sqrt(x*x + y*y)/(npix/2)
    # only select pixels above horizon
    t = (d <= 1.0)
    n_total = t.sum()
    dOMEGA.fill( math.pi*2.00/n_total )
    
    za = numpy.ones((npix,npix),dtype=numpy.float32)*numpy.pi/2.0
    
    if projection == 'SIN' :
        za[t]  =  numpy.arcsin(d[t])
        dOMEGA = numpy.cos(za)*math.pi*2.00/n_total
        print('Using slant orthographic projection')
    elif projection == 'ZEA': # https://casa.nrao.edu/aips2_docs/memos/107/node2.html#SECTION00022200000000000000
                            # https://asdf-standard.readthedocs.io/en/stable/schemas/stsci.edu/asdf/transform/zenithal_equal_area-1.0.0.html
        d     = d*2**0.5; #ZEA requires R to extend beyond 1.
        za[t] = 2*numpy.arcsin(d[t]/2.0) # = 2*arcsin( d * (sqrt(2)/2) ) = 2 * arcsin( d / sqrt(2) )
        print('Using zenithal equal area projection')
    else :
        raise Exception("Unknown projection %s" % (projection))

    az     = numpy.arctan2(y,x)
    dOMEGA_sum = dOMEGA.sum()
    print("DEBUG : dOMEGA_sum = %.8f" % (dOMEGA_sum))
    
    return (az,za,n_total,dOMEGA)


def get_fits_beam_multi( azim_rad, za_rad, frequency_mhz, 
                         polarisation="X", 
                         station_name="EDA", 
                         simulation_path="$HOME/aavs-calibration/BeamModels/",
                         use_isotropic_antenna=False,
                         use_beam_as_is=True,
                         debug=False
                        ) :

   global global_fits_counter 
   global_fits_counter = global_fits_counter + 1
#   global azim_map
#   global za_map
   
   azim_deg = azim_rad*(180.00/math.pi)
   za_deg   = za_rad*(180.00/math.pi)

   print("requestion beam model for azza map of size (%d x %d) pixels" % (azim_deg.shape[0],azim_deg.shape[1]))

   (current_fits_beam,beam_file_name) = read_beam_fits( frequency_mhz, polarisation, station_name, simulation_path )
      
   x_size = int( current_fits_beam[0].header['NAXIS1'] )
   y_size = int( current_fits_beam[0].header['NAXIS2'] )
   npix = x_size
   x_center = x_size / 2
   y_center = y_size / 2 
   print('INFO : Current fits file %s size = %d x %d -> center = %d x %d' % (beam_file_name,x_size,y_size,x_center,y_center))   
   
   # assuming that ds9 shows :   
   #       NORTH
   # EAST        WEST
   #       SOUTH 
   # FITS files are in orthographic projection :
   # d = np.sqrt(x*x + y*y)/(npix/2)
#   z = numpy.linspace(-npix/2.0,npix/2.0,num=npix) # was dtype=np.float32  
#   x = numpy.empty((npix,npix),dtype=numpy.float32)
#   y = numpy.empty((npix,npix),dtype=numpy.float32)
#   for i in range(npix):
#        y[i,0:] = z
#        x[0:,i] = z
#   d = numpy.sqrt(x*x + y*y)/(npix/2)

   # only select pixels above horizon
#   t = (d <= 1.0)
#   za = numpy.zeros((npix,npix),dtype=numpy.float32)*numpy.NaN
#   za[t] = numpy.arcsin(d[t])
#   print 'Using slant orthographic projection'
#   az = numpy.arctan2(y,x)
#   az = az+math.pi #0 to 2pi
#   az = 2*math.pi-az #Change to clockwise from top (when origin is in top-left)


# BUG IS PROBABLY HERE : !!! 
# (Pdb) p azim_deg[256,400]
# 90.19825
# (Pdb) 

   # see /home/msok/ska/aavs/aavs0.5/trunk/simulations/FEKO/beam_models/MWA_EE/MWAtools_pb/eda/beam_tools.py
   # get x,y from azim,za ... - other way around
   # HERE I AM ASSUMING THAT BeamModel fits files in ~/aavs-calibration/BeamModel/ are in SIN-orthographic projection :
   # and I am using the same mapping as passed in parameters (az,za) - this is probably wrong as one is ZAE and the other SIN - this is 
   # probably the reason for the holes !!!
   # I need a new AZZA map here for the sin projection !
   # Is it correct assumption ?
#   print "Image size = (%d,%d) vs. map size (%d,%d)" % (x_size,y_size,azim_map.shape[0],azim_map.shape[1])
#   if x_size != azim_map.shape[0] or y_size != azim_map.shape[1] :
#      raise Exception(  "ERROR : %d != %d or %d != %d -> size mismatch -> cannot continue" % (x_size,y_size,azim_map.shape[0],azim_map.shape[1]) )
      

   dist_za = None 
   if beam_file_name.upper().find("ZEA") >= 0 :
      # if files / AZZA map is in ZEA (zenith equal area projection) -> distance from the center has to by calculated in a different way 
      # possibly still division by sqrt(2.00) is needed here as otherwise will go beyond 1 
      dist_za = (npix/2)*numpy.sqrt( 2.00*(1.00 - numpy.cos(za_deg*(math.pi/180.00)) )  ) / math.sqrt(2.00)
   else :
      # just normal SIN projection :
      dist_za = numpy.sin( za_deg*(math.pi/180.00) )*(npix/2)
      
   x_pixel = ( npix/2 + dist_za*numpy.cos( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
   y_pixel = ( npix/2 + dist_za*numpy.sin( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
      
   reshape_x_size = npix*npix
   reshape_y_size = npix*npix
   if (npix*npix) != (azim_deg.shape[0]*azim_deg.shape[1] ) :
      if azim_deg.shape[0]==1 and azim_deg.shape[1]==1 :
         reshape_x_size = 1 # int( math.sqrt( azim_deg.shape[0]*azim_deg.shape[1]  ) )
         reshape_y_size = 1 # reshape_x_size
      else :
         raise Exception("(azim,za) map has to be either the same size as beam fits files (%d,%d) or (1,1) - (%d,%d) is not allowed" % (x_size,y_size,azim_deg.shape[0],azim_deg.shape[1]))
   
   print("x_pixel : count = %d, npix=%d -> reshape_x_size = %d" % (len(x_pixel),npix,reshape_x_size))
 
   beam_values = numpy.zeros( azim_deg.shape )
      
   # ISOTROPIC ANTENNA TEST :
   if use_isotropic_antenna :
      beam_values = numpy.ones( azim_deg.shape )
      return (beam_values,current_fits_beam,x_pixel.astype(int),y_pixel.astype(int))
      

   if len(azim_deg) > 1 :
       if use_beam_as_is : 
          if beam_values.shape[0] != current_fits_beam[0].data.shape[0] or beam_values.shape[1] != current_fits_beam[0].data.shape[1] :
             print("ERROR in code shapes do not match (%d != %d or %d != %d)" % (beam_values.shape[0],current_fits_beam[0].data.shape[0],beam_values.shape[1],current_fits_beam[0].data.shape[1]))
             return (None,None,None,None)
         
          beam_values = numpy.sqrt( current_fits_beam[0].data )

       else :   
           if len(azim_deg) > 1 :
              for i in range(0,len(x_pixel)) :
                 for j in range(0,len(y_pixel)):
                    if not numpy.isnan(x_pixel[i,j]) and not numpy.isnan(y_pixel[i,j]) :
                        x_int = int(x_pixel[i,j])
                        y_int = int(y_pixel[i,j])
                          
#                print "(x_int,y_int) = (%d,%d)" % (x_int,y_int)
                
                        if x_int < beam_values.shape[0] and y_int < beam_values.shape[1] :
                           beam_values[x_int,y_int] = numpy.sqrt( current_fits_beam[0].data[x_int,y_int] )
                        else :
                           print("WARNING : skipped pixel (x_int,y_int) = (%d,%d)" % (x_int,y_int))
   else :
       x_pixel_int = int(x_pixel)
       y_pixel_int = int(y_pixel)

       beam_values[0,0] = numpy.sqrt( current_fits_beam[0].data[x_pixel_int,y_pixel_int] )

#  test save :
   if debug : 
      hdu = pyfits.PrimaryHDU()
      hdu.data = beam_values
      out_fits_name = "test_beam_fits_%05d.fits" % (global_fits_counter)
      hdulist = pyfits.HDUList([hdu])
      hdulist.writeto( out_fits_name , clobber=True )
      
      hdu_az = pyfits.PrimaryHDU()
      hdu_az.data = azim_deg
      out_fits_name = "test_azim_fits_%05d.fits" % (global_fits_counter)
      hdulist_az = pyfits.HDUList([hdu_az])
      hdulist_az.writeto( out_fits_name , clobber=True )

      hdu_za = pyfits.PrimaryHDU()
      hdu_za.data = za_deg
      out_fits_name = "test_za_fits_%05d.fits" % (global_fits_counter)
      hdulist_za = pyfits.HDUList([hdu_za])
      hdulist_za.writeto( out_fits_name , clobber=True )

   

         
   
   # multiple pixels :
#   x_pixel_lin = x_pixel.reshape( reshape_x_size )
#   y_pixel_lin = y_pixel.reshape( reshape_y_size )
#   x_pixel_lin_nonan = []
#   y_pixel_lin_nonan = []
#
#   for i in range(0,len(x_pixel_lin)) :
#       if not numpy.isnan( x_pixel_lin[i] ) and not numpy.isnan( y_pixel_lin[i] ) :
#           x_pixel_lin_nonan.append( x_pixel_lin[i] )
#           y_pixel_lin_nonan.append( y_pixel_lin[i] )
#   
#   x_pixel_lin_nonan = numpy.array( x_pixel_lin_nonan )
#   y_pixel_lin_nonan = numpy.array( y_pixel_lin_nonan )
#   cnt_nonan = len( x_pixel_lin_nonan )
     
##   beam_values = current_fits_beam[0].data[y_pixel.astype(int),x_pixel.astype(int)]
#   for i in range(0,cnt_nonan) :
#      y_int = int(y_pixel_lin_nonan[i])
#      for j in range(0,cnt_nonan) :
#          x_int = int(x_pixel_lin_nonan[i])
##          print "(%d,%d) -> (%d,%d) = %.4f" % (i,j,x_int,y_int,current_fits_beam[0].data[y_int,x_int])
#          beam_values[x_int,y_int] = current_fits_beam[0].data[x_int,y_int] # or x,y ???

#   if len(x_pixel_lin_nonan) != 1 and len(y_pixel_lin_nonan) != 1 :
#      beam_values[x_pixel_lin_nonan.astype(int),y_pixel_lin_nonan.astype(int)] = current_fits_beam[0].data[x_pixel_lin_nonan.astype(int),y_pixel_lin_nonan.astype(int)] 
#   else :
#      beam_values[0,0] = current_fits_beam[0].data[x_pixel_lin_nonan.astype(int),y_pixel_lin_nonan.astype(int)]
      
   return (beam_values,current_fits_beam,x_pixel.astype(int),y_pixel.astype(int))


if __name__ == "__main__":

   do_remapping = False
   
   if do_remapping :
      fits_name = "EDA_Xpol_ortho_204.fits"
      if len(sys.argv) > 1:
         fits_name = sys.argv[1]
   
      step = 1 
      if len(sys.argv) > 2:
         step = int( sys.argv[2] )
  
      
      print("Remapping fits_file = %s in steps of %d" % (fits_name,step))
      remap_beam( fits_name , step=step, do_test=False, radius=20 )
   else :
      az = 0 
      za = 0
      freq_mhz = 160
      
      if len(sys.argv) > 1:
         az = float( sys.argv[1] )

      if len(sys.argv) > 2:
         za = float( sys.argv[2] )

      if len(sys.argv) > 3:
         freq_mhz = float( sys.argv[3] )

      beam_x = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='X' )   
      beam_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='Y' )   
      
      print("BEAM_X = %.4f , BEAM_Y = %.4f " % (beam_x,beam_y))
      