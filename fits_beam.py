from __future__ import print_function
# from pylab import *
# import pdb

import numpy
import os,sys
import math
from string import Template
import datetime
import time
import beam_tools
import copy

# in order to avoid crash due to missing SUPER-EXACT ephemeris files which will resule in unnoticeable error:
import astropy
from astropy.utils.iers import conf
conf.auto_max_age = None
from astropy.utils import iers
iers.conf.auto_download = False  


# import pyfits
import astropy.io.fits as pyfits
try :
   from astropy.coordinates import SkyCoord, EarthLocation, get_sun
   # CONSTANTS :
   MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)
   
   import sky2pix
   from astropy.time import Time
except :
   print("ERROR : could not load astropy.coordinates - fits2beam cannot be used")   
   MWA_POS=None

try :
   import h5py 
except :
   print("WARNING : could not load package h5py -> option --outfile_beam_on_sun will not be available, but all the others are ok")

# options :
from optparse import OptionParser,OptionGroup

# current_fits_filename = None
current_fits_beams    = {}
# azim_map              = None
# za_map                = None


def read_list( list_filename ) :
   out_fits_list=[]

   print("DEBUG : read_list : reading file %s" % (list_filename))

   if os.path.exists( list_filename ) and os.stat( list_filename ).st_size > 0 :
      file=open( list_filename ,'r')
      data=file.readlines()
      for line in data :
         line=line.rstrip()
         words = line.split(' ')
         if line[0] == '#' :
            continue

         if line[0] != "#" :
            elem = words[0+0] 
            out_fits_list.append( elem )
                        
      file.close()
   else :
      print("WARNING : empty or non-existing file %s" % (list_filename))

   return out_fits_list



def save_fits( data , out_fits_name ) :
   hdu = pyfits.PrimaryHDU()
   hdu.data = data
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_fits_name , overwrite=True )

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
         
def matlab2sin( matlab_fits_file , x_size=512, step=1, do_test=False, radius=20 ) :
   beam   = pyfits.open( matlab_fits_file )
   matlab_data = beam[1] # was 0 previously
   y_size = x_size
   npix = x_size    
#   print("Read file %s with size (%d,%d)" % (matlab_fits_file,x_size,y_size))
   if matlab_data is None :
       print("ERROR : could not read FITS file %s" % (matlab_fits_file))
       return 

   print("Read file %s with size (%d,%d)" % (matlab_fits_file,matlab_data.shape[0],matlab_data.shape[0]))
   
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
   
   out_file = matlab_fits_file.replace(".fits","_test.fits" )
   save_fits( beam[0].data, out_file )
      
   out_beam = numpy.zeros( (x_size,y_size) )
   for y in range(0,y_size,step):
      # print "i = %d" % (i)
      for x in range(0,x_size): 
         print("TEST (x,y) = (%d,%d)" % (x,y))
         az_zea_ij = az_zea[y,x] # this order of x,y is ok (see makeAZZA)
         za_zea_ij = za_zea[y,x] # this order of x,y is ok (see makeAZZA)
         
         if not numpy.isnan( za_zea_ij ) :
             d_sin     = math.sin( za_zea_ij )
             y_zea_ij  = y_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
             x_zea_ij  = x_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
         
             d_sin = (npix/2)*math.sin( za_zea_ij )
             y_sin_ij = d_sin*math.cos( az_zea_ij ) + (npix/2)
             x_sin_ij = d_sin*math.sin( az_zea_ij ) + (npix/2)
             
             if do_test :     
                 print("\t(%d,%d) : ZEA (x_zea,y_zea) = (%d,%d) , (az,za) = (%.2f,%.2f) [deg] -> SIN : (x_sin,y_sin) = (%d,%d)" % (x,y,x_zea_ij,y_zea_ij,az_zea_ij,za_zea_ij,int(x_sin_ij),int(y_sin_ij)))

             # find given (az,za) from ZEA map in SIN map and return beam value (in SIN mapping) :
             beam_value = numpy.NaN
             
             az_zea_ij_deg = az_zea_ij*(180.00/math.pi)
             za_zea_ij_deg = za_zea_ij*(180.00/math.pi)
             
             az_zea_ij_deg_idx = int( az_zea_ij_deg / 0.5 )
             za_zea_ij_deg_idx = int( za_zea_ij_deg / 0.5 )
             
             beam_value = matlab_data.data[az_zea_ij_deg_idx,za_zea_ij_deg_idx]
             out_beam[x,y] = beam_value
                  

   out_file = matlab_fits_file.replace(".fits","_aee.fits" )    
   save_fits( out_beam, out_file )            
         
def beammap2sin( theta_phi, beam_values_2d , x_size=512, step=1, do_test=False, radius=20, out_file="feko.fits" ) :
   y_size = x_size
   npix = x_size

   print("Read file with size (%d,%d)" % (theta_phi.shape[0],theta_phi.shape[1]))
   
   (az_sin,za_sin,x_sin,y_sin,z_sin,d_sin) = beam_tools.makeAZZA( x_size, 'SIN' ,azim_from_north=True, return_all=True, force_zenith=True )
   (az_zea,za_zea,x_zea,y_zea,z_zea,d_zea) = beam_tools.makeAZZA( x_size, 'ZEA' ,azim_from_north=True, return_all=True, force_zenith=True )
   
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
   
   test_out_file = out_file.replace(".fits","_test.fits" )
   save_fits( beam_values_2d, test_out_file )
   
   out_beam = numpy.zeros( (x_size,y_size) )
   for y in range(0,y_size-1,step):
      # print "i = %d" % (i)
      for x in range(0,x_size-1): 
         # print("TEST (x,y) = (%d,%d)" % (x,y))
         az_zea_ij = az_zea[y,x] # this order of x,y is ok (see makeAZZA)
         za_zea_ij = za_zea[y,x] # this order of x,y is ok (see makeAZZA)
         
         if not numpy.isnan( za_zea_ij ) :
             d_sin     = math.sin( za_zea_ij )
             y_zea_ij  = y_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
             x_zea_ij  = x_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
         
             d_sin = (npix/2)*math.sin( za_zea_ij )
             y_sin_ij = d_sin*math.cos( az_zea_ij ) + (npix/2)
             x_sin_ij = d_sin*math.sin( az_zea_ij ) + (npix/2)
             
             if do_test :     
                 print("\t(%d,%d) : ZEA (x_zea,y_zea) = (%d,%d) , (az,za) = (%.2f,%.2f) [deg] -> SIN : (x_sin,y_sin) = (%d,%d)" % (x,y,x_zea_ij,y_zea_ij,az_zea_ij,za_zea_ij,int(x_sin_ij),int(y_sin_ij)))

             # find given (az,za) from ZEA map in SIN map and return beam value (in SIN mapping) :
             beam_value = numpy.NaN
             
             az_zea_ij_deg = az_zea_ij*(180.00/math.pi)
             za_zea_ij_deg = za_zea_ij*(180.00/math.pi)
             
             az_zea_ij_deg_idx = int( round(az_zea_ij_deg / 0.5) )
             za_zea_ij_deg_idx = int( round(za_zea_ij_deg / 0.5) )
             
             if az_zea_ij_deg_idx < beam_values_2d.shape[0] and za_zea_ij_deg_idx < beam_values_2d.shape[1] : 
                beam_value = beam_values_2d[az_zea_ij_deg_idx,za_zea_ij_deg_idx]
#                if (x+1) < out_beam.shape[0] and (y+1) < out_beam.shape[1] :
                out_beam[x,y] = beam_value
                
                if x>=250 and x<=260 and y>=250 and y<=260 :
                   print("DEBUG(%d,%d) = %.8f at (az,za) = (%.6f,%.6f)" % (x,y,beam_value,az_zea_ij,za_zea_ij))
                  

   save_fits( out_beam, out_file )            
         
def beammap2sin_OLD( theta_phi, beam_values_2d , x_size=512, step=1, do_test=False, radius=20, out_file="feko.fits" ) :
   y_size = x_size
   npix = x_size

   print("Read file with size (%d,%d)" % (theta_phi.shape[0],theta_phi.shape[1]))
   
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
   
   test_out_file = out_file.replace(".fits","_test.fits" )
   save_fits( beam_values_2d, test_out_file )
   
   out_beam = numpy.zeros( (x_size,y_size) )
   for y in range(0,y_size-1,step):
      # print "i = %d" % (i)
      for x in range(0,x_size-1): 
         # print("TEST (x,y) = (%d,%d)" % (x,y))
         az_zea_ij = az_zea[y,x] # this order of x,y is ok (see makeAZZA)
         za_zea_ij = za_zea[y,x] # this order of x,y is ok (see makeAZZA)
         
         if not numpy.isnan( za_zea_ij ) :
             d_sin     = math.sin( za_zea_ij )
             y_zea_ij  = y_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
             x_zea_ij  = x_zea[x,y] + (npix/2) # or [y,x] see "Test orientation xy SIN :"
         
             d_sin = (npix/2)*math.sin( za_zea_ij )
             y_sin_ij = d_sin*math.cos( az_zea_ij ) + (npix/2)
             x_sin_ij = d_sin*math.sin( az_zea_ij ) + (npix/2)
             
             if do_test :     
                 print("\t(%d,%d) : ZEA (x_zea,y_zea) = (%d,%d) , (az,za) = (%.2f,%.2f) [deg] -> SIN : (x_sin,y_sin) = (%d,%d)" % (x,y,x_zea_ij,y_zea_ij,az_zea_ij,za_zea_ij,int(x_sin_ij),int(y_sin_ij)))

             # find given (az,za) from ZEA map in SIN map and return beam value (in SIN mapping) :
             beam_value = numpy.NaN
             
             az_zea_ij_deg = az_zea_ij*(180.00/math.pi)
             za_zea_ij_deg = za_zea_ij*(180.00/math.pi)
             
             az_zea_ij_deg_idx = int( round(az_zea_ij_deg / 0.5) )
             za_zea_ij_deg_idx = int( round(za_zea_ij_deg / 0.5) )
             
             if az_zea_ij_deg_idx < beam_values_2d.shape[0] and za_zea_ij_deg_idx < beam_values_2d.shape[1] : 
                beam_value = beam_values_2d[az_zea_ij_deg_idx,za_zea_ij_deg_idx]
#                if (x+1) < out_beam.shape[0] and (y+1) < out_beam.shape[1] :
                out_beam[x,y] = beam_value
                  

   save_fits( out_beam, out_file )            
         
         


# TODO (2020-08-20) : use dictionary to cache all FITS files ! It all works very slow on bighorns ...
def read_beam_fits( frequency_mhz, polarisation="X", station_name="EDA", simulation_path="$HOME/aavs-calibration/BeamModels/" , postfix="zea" ) :
#   global current_fits_filename
   global current_fits_beams
#   global azim_map
#   global za_map

   polarisation_string = polarisation[0]
   postfix_full = ""
   if postfix is not None and len(postfix)>0 :
      postfix_full = "_" + postfix 
   beam_file_name_template = "%s/%s/%s/%s_%spol_ortho_%03d%s.fits" % (simulation_path,station_name,postfix.upper(),station_name,polarisation_string,int(frequency_mhz),postfix_full)
   beam_file_name = Template( beam_file_name_template ).substitute(os.environ)
   fits_beam = None
   
#   if current_fits_filename is None or current_fits_beam is None or beam_file_name != current_fits_filename :
   cached=False
   if beam_file_name in current_fits_beams.keys() :
      print("FITS file %s found in cache -> returning" % (beam_file_name))
      cached=True
   else :
      print("Reading fits file %s ..." % (beam_file_name))
           
      fits_beam = pyfits.open( beam_file_name )
      x_size = fits_beam[0].header['NAXIS1']  
      y_size = fits_beam[0].header['NAXIS2']

      print() 
      print('# \tRead fits file %s of size %d x %d' % (beam_file_name,x_size,y_size))
      
      if len(current_fits_beams) < 50 :
         current_fits_beams[beam_file_name] = copy.copy( fits_beam )
         cached=True
      else :
         print("INFO : exceeded number of cached beams -> no caching for file %s" % (beam_file_name))
#      fits_beam.close()
#      current_fits_filename = beam_file_name
      
      
#   x_size = current_fits_beam[0].header['NAXIS1']  
#   y_size = current_fits_beam[0].header['NAXIS2']
#   if azim_map is None or za_map is None :
#      (azim_map,za_map) = makeAZZA( x_size )
#      print "Generated local AZZA map for size %d x %d pixels" % (x_size,x_size)
  
   if cached :
      return (current_fits_beams[beam_file_name],beam_file_name)
   else :
      return (fits_beam,beam_file_name)

# individual_antenna_beams_path="~/aavs-calibration/
# EDA_Xpol_ortho_169.fits
def get_fits_beam( azim_deg, za_deg, frequency_mhz, polarisation="X", station_name="EDA", simulation_path="$HOME/aavs-calibration/BeamModels/", projection="zea", debug_level=0 ) :
   print("get_fits_beam : requestion beam model for azza map of size (%d x %d) pixels , projection = %s" % (azim_deg.shape[0],azim_deg.shape[1],projection))

   (current_fits_beam,beam_file_name) = read_beam_fits( frequency_mhz, polarisation, station_name, simulation_path, postfix=projection )
   

   x_size = int( current_fits_beam[0].header['NAXIS1'] )
   y_size = int( current_fits_beam[0].header['NAXIS2'] )
   npix = x_size
   x_center = x_size / 2
   y_center = y_size / 2 
   if debug_level > 0 :
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
   
   
#   dist_za = numpy.sin( za_deg*(math.pi/180.00) )*(npix/2)
# 2020-05-30 fix :   
   is_zea = False
   if beam_file_name.upper().find("ZEA") >= 0 or beam_file_name.upper().find("AEE") >= 0 : # AEE also in ZEA projection !   
      is_zea = True

   dist_za = None 
   if is_zea :
      # if files / AZZA map is in ZEA (zenith equal area projection) -> distance from the center has to by calculated in a different way 
      # possibly still division by sqrt(2.00) is needed here as otherwise will go beyond 1 
      dist_za = (npix/2)*numpy.sqrt( 2.00*(1.00 - numpy.cos(za_deg*(math.pi/180.00)) )  ) / math.sqrt(2.00)
   else :
      # just normal SIN projection :
      dist_za = numpy.sin( za_deg*(math.pi/180.00) )*(npix/2)

   
   x_pixel = ( npix/2 + dist_za*numpy.cos( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
   y_pixel = ( npix/2 + dist_za*numpy.sin( azim_deg*(math.pi/180.00) ) )  # orientaition in python is N(az=0) = (512,256), E(az=90deg) = (256,512) , S(az=180deg) = (0,255), W(az=-90/270deg) = (255,0)
   if debug_level > 0 :
      print("x_pixel : count = %d" % (len(x_pixel)))

   # single pixel :
#   beam_value = current_fits_beam[0].data[int(y_pixel),int(x_pixel)]
   x_int = int( round(x_pixel) )
   y_int = int( round(y_pixel) )
   
   if x_int >= 0 and x_int < x_size and y_int >=0 and y_int < y_size : 
      try : 
         beam_value = current_fits_beam[0].data[ x_int , y_int ] # 2020-05-30 : changed to the same orientation as in get_fits_beam_multi - uses the fact of transposition in python so that the answer is correct !
                                                                 #              otherwise it seems to have been giving wrong answers !!!
                                                                 #              no x <-> y swap here to return transposed value (North is on top in ds9) !
      except :
         print("ERROR : exception caught when accessing data at pixel (%d,%d) corresponding to (azim,za) = (%.4f,%.4f) [deg] -> skipped" % (x_int , y_int, azim_deg, za_deg))
   else :      
      print("WARNING : calculated pixel coordinates for (azim,za) = (%.4f,%.4f) [deg] are (%d,%d) which is outside the image -> asigning closest value" % (azim_deg, za_deg, x_int , y_int))
      
      if x_int <=0 :
         x_int = 0
      if y_int <= 0 :
         y_int = 0
      if x_int >= x_size :
         x_int = (x_size-1)
      if y_int >= y_size :
         y_int = (y_size-1)
         
      print("WARNING : getting beam value at pixel (%d,%d)" % (x_int,y_int))   
         
      beam_value = current_fits_beam[0].data[ x_int , y_int ] # 2020-05-30 : changed to the same orientation as in get_fits_beam_multi - uses the fact of transposition in python so that the answer is correct 
                                                              #              otherwise it seems to have been giving wrong answers !!!
                                                              #              no x <-> y swap here to return transposed value (North is on top in ds9) !
                                                                                       
   print("Beam( %.4f deg, %.4f deg ) = %.8f at pixel (%d,%d) , dist_za = %.4f [deg]" % (azim_deg,za_deg,beam_value,x_int,y_int,dist_za))
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
                         debug=False, 
                         power=False,
                         projection="zea"
                        ) :
   print("DEBUG get_fits_beam_multi : projection = %s" % (projection))

   global global_fits_counter 
   global_fits_counter = global_fits_counter + 1
#   global azim_map
#   global za_map
   
   azim_deg = azim_rad*(180.00/math.pi)
   za_deg   = za_rad*(180.00/math.pi)

   print("get_fits_beam_multi : requestion beam model for azza map of size (%d x %d) pixels, projection = %s" % (azim_deg.shape[0],azim_deg.shape[1],projection))

   (current_fits_beam,beam_file_name) = read_beam_fits( frequency_mhz, polarisation, station_name, simulation_path, postfix=projection )
      
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
 
   is_zea = False
   if beam_file_name.upper().find("ZEA") >= 0 or beam_file_name.upper().find("AEE") >= 0 : # AEE also in ZEA projection !   
      is_zea = True            

   dist_za = None 
   if is_zea :
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

          beam_values = None 
          if power :
             beam_values = current_fits_beam[0].data
          else :         
             beam_values = numpy.sqrt( current_fits_beam[0].data )

       else :   
           if len(azim_deg) > 1 :
              for i in range(0,len(x_pixel)) :
                 for j in range(0,len(y_pixel)):
                    if not numpy.isnan(x_pixel[i,j]) and not numpy.isnan(y_pixel[i,j]) :
                        x_int = int( round(x_pixel[i,j]) )
                        y_int = int( round(y_pixel[i,j]) )
                          
#                print "(x_int,y_int) = (%d,%d)" % (x_int,y_int)
                
                        if x_int < beam_values.shape[0] and y_int < beam_values.shape[1] :
                           if power :
                              # power is no SQRT :
                              beam_values[x_int,y_int] = current_fits_beam[0].data[x_int,y_int]
                           else :
                              # Jones is SQRT :
                              beam_values[x_int,y_int] = numpy.sqrt( current_fits_beam[0].data[x_int,y_int] ) # same orientation as in get_fits_beam - uses the fact of transposition in python so that the answer is correct 
                                                                                                              # no x <-> y swap here to return transposed value (North is on top in ds9) !
                        else :
                           print("WARNING : skipped pixel (x_int,y_int) = (%d,%d)" % (x_int,y_int))
   else :
       x_pixel_int = int( round(x_pixel) )
       y_pixel_int = int( round(y_pixel) )

       if power :
          beam_values[0,0] = current_fits_beam[0].data[x_pixel_int,y_pixel_int]
       else : 
          beam_values[0,0] = numpy.sqrt( current_fits_beam[0].data[x_pixel_int,y_pixel_int] ) # same orientation as in get_fits_beam - uses the fact of transposition in python so that the answer is correct 
                                                                                           # no x <-> y swap here to return transposed value (North is on top in ds9) !

#  test save :
   if debug or False : 
      hdu = pyfits.PrimaryHDU()
      hdu.data = beam_values
      out_fits_name = "test_beam_fits_%05d.fits" % (global_fits_counter)
      hdulist = pyfits.HDUList([hdu])
      hdulist.writeto( out_fits_name , overwrite=True )
      
      hdu_az = pyfits.PrimaryHDU()
      hdu_az.data = azim_deg
      out_fits_name = "test_azim_fits_%05d.fits" % (global_fits_counter)
      hdulist_az = pyfits.HDUList([hdu_az])
      hdulist_az.writeto( out_fits_name , overwrite=True )

      hdu_za = pyfits.PrimaryHDU()
      hdu_za.data = za_deg
      out_fits_name = "test_za_fits_%05d.fits" % (global_fits_counter)
      hdulist_za = pyfits.HDUList([hdu_za])
      hdulist_za.writeto( out_fits_name , overwrite=True )

   

         
   
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

def beam_correct_flux( dt_arr, tm_arr, az_arr, el_arr, ra_arr, dec_arr, flux_arr, cnt, frequency_mhz, outfile="out.txt", pol="X", projection="zea" ) :
   print("beam_correct_flux :")
   
   out_f = open( outfile , "w" )
   line = ("# DATE TIME AZ[deg] EL[deg] RA[deg] DEC[deg] Flux[Jy] Beam_%s FluxCorr[Jy]\n" % (pol))
   out_f.write( line )
   for i in range(0,cnt) :
      dt = dt_arr[i]
      tm = tm_arr[i]
      az = az_arr[i]
      el = el_arr[i]
      za = 90.00 - el
      ra = ra_arr[i]
      dec = dec_arr[i]
      flux = flux_arr[i]

      beam_pol = 1.00      
      if pol == "I" :
         # this is not correct way as X and Y should be beam-corrected before averaging !
         
         beam_x = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , frequency_mhz=frequency_mhz, polarisation="X", projection=projection )
         beam_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , frequency_mhz=frequency_mhz, polarisation="Y", projection=projection )
         beam_pol = (beam_x + beam_y)/2.00
         
         print("WARNING : pol=%s will use Beam_I = (X+Y)/2 = (%.8f + %.8f)/2.00 = %.8f, but X and Y images should be beam-corrected separately first !" % (pol,beam_x,beam_y,beam_pol))
      else :
         beam_pol = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , frequency_mhz=frequency_mhz, polarisation=pol, projection=projection )
      
      
      
      flux_corr = flux / beam_pol
      
      line = "%s %s %.14f %.14f %.14f %.14f %.14f %.4f %.14f\n" % (dt,tm,az,el,ra,dec,flux,beam_pol,flux_corr)
      out_f.write( line )

      
   out_f.close()

def read_time_azh_file( filename, options,
#                        dt_index   = 0, 
#                        tm_index   = 1,
#                        az_index   = 2,
#                        el_index   = 3,
#                        ra_index   = 4,
#                        dec_index  = 5,
#                        flux_index = 6,
                        
                        default_ra  = 0.00,
                        default_dec = 0.00,
                        
                        debug_level = 1
                      ) :
# format : date; time; azimuth; elevation; right ascension; declination; and peak flux density
   print("read_data(%s) ..." % (filename))
   file=open(filename,'r')

   # reads the entire file into a list of strings variable data :
   data=file.readlines()
   # print data

   if debug_level > 0 :
      print("DEBUG : %d / %d / %d / %d / %d / %d / %d" % (options.dt_index,options.tm_index,options.az_index,options.el_index,options.ra_index,options.dec_index,options.flux_index))

   # initialisation of empty lists :
   dt=[]
   tm=[]
   az=[]
   el=[]
   ra=[]
   dec=[]
   flux=[]
   cnt=0
   for line in data : 
      words = line.split(' ')

      if line[0] == '#' or line[0]=='\n' or len(line) <= 0 or len(words)<4 :
         continue
      
#      print "line = %s , words = %d" % (line,len(words))

      if line[0] != "#" :
#         print line[:-1]
# format : 2020-01-31 04:26:28.9 152.32707214144972 42.232977712074245 1.26869122140846 -62.705127298233535 142.89981842041016     
         dt1=words[options.dt_index+0]
         tm1=words[options.tm_index+0]
         az1=float(words[options.az_index+0])         
         el1=float(words[options.el_index+0])
         
         ra1 = default_ra
         if options.ra_index >= 0 :
            ra1 = float(words[options.ra_index+0])
         
         dec1 = default_dec
         if options.dec_index >= 0 :      
            dec1 = float(words[options.dec_index+0])
            
         f = float(words[options.flux_index+0])

         # format : date; time; azimuth; elevation; right ascension; declination; and peak flux density
         dt.append(dt1)
         tm.append(tm1)
         az.append(az1)
         el.append(el1)
         ra.append(ra1)
         dec.append(dec1)
         flux.append(f)
         
         if debug_level > 0 :
            print("DEBUG : %s / %s (az,el) = (%.4f,%.4f) [deg], (ra,dec) = (%.4f,%.4f) [deg] , flux = %.2f [Jy]" % (dt1,tm1,az1,el1,ra1,dec1,f))
         
         cnt += 1
         

   print("Read %d lines from file %s" % (cnt,filename))

   return (dt,tm,az,el,ra,dec,flux,cnt)


def save_beam_on_sun_file( options ) :
   out_file = "beam_on_sun.txt"
   if options.outfile_beam_on_sun is not None :
      out_file = options.outfile_beam_on_sun
      
   if options.infile_hdf5list is None or not os.path.exists( options.infile_hdf5list ) :
      print("ERROR : input file with list of hdf5 files not specified (None) or %s does not exist -> cannot continue" % (options.infile_hdf5list) )
      os.sys.exit(-1)

   (hdf5_files_list) = read_list( options.infile_hdf5list )
   
   print("DEBUG : read %d files from list file %s" % (len(hdf5_files_list),options.infile_hdf5list))
   
   if len(hdf5_files_list) > 0 :
   
      out_f = open( options.outfile_beam_on_sun , "w" )
      out_f.write( "# Channel SUN_BEAM_X SUN_BEAM_Y UXTIME FREQ_MHz HDF5-FILE\n" )
      
      for hdf5_file in hdf5_files_list :
         hdf5_f = h5py.File( hdf5_file )         
         ux = hdf5_f['sample_timestamps']['data'][0]
         channel_id = hdf5_f['root'].attrs['channel_id']
         freq_mhz = float( channel_id )*(400.00/512.00)
         
         if freq_mhz >= 49 and freq_mhz <= 350 : # only these FITS files are available for the beam models 
            ( ra, dec, az, alt, za ) = sun_position( ux )
            print("INFO : %s -> ux = %.4f , channel_id = %d (%.2f MHz) -> sun position calculated (RA,DEC) = (%.4f,%.4f) [deg] , (AZ,ELEV,ZA) = (%.4f,%.4f,%.4f) [deg]" % (hdf5_file, ux, channel_id, freq_mhz, ra, dec, az, alt, za))
         
            beam_x = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='X', projection=options.projection, station_name=options.station_name )   
            beam_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='Y', projection=options.projection, station_name=options.station_name )   
            print("\tBEAM_X = %.4f , BEAM_Y = %.4f " % (beam_x,beam_y))
         
            line = "%d %.4f %.4f %.4f %.4f %s\n" % (channel_id,beam_x,beam_y,ux,freq_mhz,hdf5_file)
            out_f.write( line )
         else :
            print("WARNING : channel_id = %d , freq = %.2f MHz - skipped (missing beam model FITS files)" % (channel_id,freq_mhz))

      out_f.close()   
         
   else :
      print("ERROR : no files in the list file %s" % (hdf5_files_list))
      os.sys.exit(-1)

def fits2beam( fitsfile, options=None ) :

   print("Reading fits file %s" % (fitsfile))
   fits = pyfits.open(fitsfile)
   x_size=fits[0].header['NAXIS1']
   y_size=fits[0].header['NAXIS2']
   dateobs=fits[0].header['DATE-OBS']

   ( ra , dec ) = sky2pix.pix2sky( fits, fitsfile )


   out_ra_fits = fitsfile.replace(".fits","_ra.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = ra.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_ra_fits ,overwrite=True)

   out_dec_fits = fitsfile.replace(".fits","_dec.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = dec.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_dec_fits ,overwrite=True)
      
   coord = SkyCoord( ra, dec, equinox='J2000',frame='icrs', unit='deg')
   coord.location = MWA_POS
   utc=fitsfile[9:24]  
   uxtime = time.mktime(datetime.datetime.strptime(utc, "%Y%m%dT%H%M%S").timetuple()) + 8*3600 # just for Perth !!!
   coord.obstime = Time( uxtime, scale='utc', format="unix" )
   altaz = coord.transform_to('altaz')
   az, alt = altaz.az.deg, altaz.alt.deg
   za = 90.00 - alt

   out_az_fits = fitsfile.replace(".fits","_az.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = az.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_az_fits ,overwrite=True)

   out_alt_fits = fitsfile.replace(".fits","_alt.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = alt.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_alt_fits ,overwrite=True)

   out_za_fits = fitsfile.replace(".fits","_za.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = za.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_za_fits ,overwrite=True)

#   beam_value_x = get_fits_beam_multi( az , za , options.frequency_mhz, polarisation='X', station_name=options.station_name ) # projection=options.projection, station_name=options.station_name )
#   beam_value_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='Y', projection=options.projection, station_name=options.station_name )

   beam_value_x = numpy.ones( (x_size,y_size) )*numpy.nan
   beam_value_y = numpy.ones( (x_size,y_size) )*numpy.nan
   
   for y in range(0,y_size) : 
      for x in range(0,x_size) :
         az_deg = az[y,x]
         za_deg = za[y,x]
         el_deg = alt[y,x]

         if el_deg >= 0 :
            if not numpy.isnan(az_deg) and not numpy.isnan(za_deg) :
               beam_x = get_fits_beam( numpy.array([[az_deg]]) , numpy.array([[za_deg]]) , options.frequency_mhz, polarisation='X', projection=options.projection, station_name=options.station_name )   
               beam_y = get_fits_beam( numpy.array([[az_deg]]) , numpy.array([[za_deg]]) , options.frequency_mhz, polarisation='Y', projection=options.projection, station_name=options.station_name )
               
               beam_value_x[y,x] = beam_x
               beam_value_y[y,x] = beam_y
               
#         beam_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='Y', projection=options.projection, station_name=options.station_name )   

      
#         beam_value_x[y,x] = beam_x

   out_bx_fits = fitsfile.replace(".fits","_BeamX.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = beam_value_x.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_bx_fits , overwrite=True)

   out_by_fits = fitsfile.replace(".fits","_BeamY.fits")
   hdu = pyfits.PrimaryHDU()
   hdu.data = beam_value_y.transpose()
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto( out_by_fits , overwrite=True)
   
   print("BEAM FILES written to files : %s and %s" % (out_bx_fits,out_by_fits))

def sun_position( unix_time=None ):
   if unix_time is None or unix_time <= 0 :
      print("ERROR : unix_time parameter not provided to sun_position function -> cannot continue")
      return (None,None,None,None)       
      
   t = Time(unix_time,format='unix')
   sunpos = get_sun(t)
   print("INFO : sun position = (RA,DEC) = (%.4f,%.4f) [deg]" % (sunpos.ra.degree, sunpos.dec.degree))
   
   # calculate az,za of the sun :
   coord = SkyCoord( sunpos.ra.degree, sunpos.dec.degree, equinox='J2000',frame='icrs', unit='deg')
   coord.location = MWA_POS
   coord.obstime = Time( unix_time, scale='utc', format="unix" )
   altaz = coord.transform_to('altaz')
   az, alt = altaz.az.deg, altaz.alt.deg
   za = 90.00 - alt
   
   return ( sunpos.ra.degree, sunpos.dec.degree, az, alt, za )
   


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tAccesses SKA-Low dipole files to get beam values\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-s','--station','--station_name',dest="station_name",default="EDA",help="Station name [EDA2 or AAVS2 or AAVS1] [default %default]",metavar="STRING")
   parser.add_option('-r','--remap','--do_remapping','--remapping',dest="do_remapping",default=False,action="store_true", help="Re-mapping [default %default]",metavar="STRING")
   parser.add_option('--matlab2fits',dest="matlab2fits",default=None,help="Convert rectangular MATLAB file into SIN or ZEA projection")
   parser.add_option('-p','--pol','--polarisation',dest="polarisation",default=None, help="Polarisation [default %default]")
   parser.add_option('-a','--azim','--az','--az_deg','--azim_deg',dest="azim_deg",default=0, help="Azimuth [deg]",type="float")
   parser.add_option('-z','--zenith_angle','--za','--za_deg',dest="za_deg",default=0, help="Zenith angle [deg]",type="float")
   parser.add_option('-e','--elev_deg','--el','--el_deg',dest="el_deg",default=None, help="Elevation [deg]",type="float")
   parser.add_option('-f','--freq_mhz','--frequency_mhz','--frequency',dest="frequency_mhz",default=160, help="Frequency [MHz]",type="float")
   parser.add_option('--projection',dest="projection",default="zea", help="Projection [default %default]")
   parser.add_option('--time_azh_file',dest="time_azh_file",default=None, help="File name to convert (AZ,H) [deg] -> Beam X/Y values and multiply flux if available [format :  date; time; azimuth; elevation; right ascension; declination; and peak flux density]")   
   parser.add_option('--fits2beam', dest="fits2beam", default=None, help="Generates beam for a specified FITS file [default %default]",metavar="STRING")
   parser.add_option('--ux','--uxtime','--unix_time','--unixtime',dest="unix_time",default=None, help="Unix time",type="float")
   
   # sun in the beam :
   parser.add_option('--sun','--beam_on_sun','--sun_beam',dest="beam_on_sun",default=False,action="store_true", help="Calculate beam value on the Sun [default %default]")
   parser.add_option('--outfile_beam_on_sun',dest="outfile_beam_on_sun",default=None, help="Beam on Sun output text file in format : Channel SUN_BEAM_X SUN_BEAM_Y UXTIME FREQ_MHz HDF5-FILE [default %default]")
   parser.add_option('--infile_hdf5list',dest="infile_hdf5list",default=None, help="Input list of HDF5 files to calculate value of the beam in the direction of the Sun [default %default]")
   
   # reading text file for beam correction :
   parser.add_option('--lightcurve_file' , dest="lightcurve_file",default=None, help="Lightcurve text file output from dump_pixel_radec.py")
   parser.add_option("--dt_index"        , dest="dt_index"   , default=0, help="DT index"  , type="int")
   parser.add_option("--tm_index"        , dest="tm_index"   , default=1, help="TM index"  , type="int")
   parser.add_option("--az_index"        , dest="az_index"   , default=2, help="AZ index"  , type="int")
   parser.add_option("--el_index"        , dest="el_index"   , default=3, help="EL index"  , type="int")
   parser.add_option("--ra_index"        , dest="ra_index"   , default=4, help="RA index"  , type="int")
   parser.add_option("--dec_index"       , dest="dec_index"  , default=5, help="DEC index" , type="int")
   parser.add_option("--flux_index"      , dest="flux_index" , default=6, help="Flux index", type="int")
   parser.add_option('--ra','--ra_deg'   , dest="ra_deg",default=0, help="RA [deg] - for lightcurve of RA,DEC object",type="float")
   parser.add_option('--dec','--dec_deg' , dest="dec_deg",default=0, help="DEC [deg] - for lightcurve of RA,DEC object",type="float")
      
   (options, args) = parser.parse_args(sys.argv[idx:])
   
   if options.station_name.upper() == "AAVS2" :
      options.station_name = "SKALA4"

   if options.station_name.upper() == "AAVS1" :
      options.station_name = "SKALA2"

   return (options, args)



if __name__ == "__main__":

   (options, args) = parse_options()
   
   if options.el_deg is not None :
      print("DEBUG : elevation parameter used = %.4f [deg]" % (options.el_deg))
      options.za_deg = (90.00 - options.el_deg)
   
   print("######################################################")
   print("PARAMETERS :")
   print("######################################################")
   print("Projection      = %s" % (options.projection))
   print("lightcurve_file = %s" % (options.lightcurve_file))
   print("time_azh_file   = %s" % (options.time_azh_file))
   print("Coordinates (az,za) = (%.4f,%.4f) [deg] ( elevation param = %s )" % (options.azim_deg,options.za_deg,options.el_deg))
   print("fits2beam       = %s" % (options.fits2beam))
   print("matlab2fits     = %s" % (options.matlab2fits))
   print("Beam on Sun values:")
   print("\tBeam on Sun         = %s" % (options.beam_on_sun))
   print("\tinfile_hdf5list     = %s" % (options.infile_hdf5list))
   print("\toutfile_beam_on_sun = %s" % (options.outfile_beam_on_sun))
   print("######################################################")

   if options.do_remapping :
      fits_name = "EDA_Xpol_ortho_204.fits"
      if len(sys.argv) > 1:
         fits_name = sys.argv[1]
   
      step = 1 
      if len(sys.argv) > 2:
         step = int( sys.argv[2] )
  
      
      print("Remapping fits_file = %s in steps of %d" % (fits_name,step))
      remap_beam( fits_name , step=step, do_test=False, radius=20 )

   elif options.time_azh_file is not None :
      print("Reading file %s ..." % (options.time_azh_file))
      (dt,tm,az,el,ra,dec,flux,cnt) = read_time_azh_file( options.time_azh_file, options ) 
      projection = options.projection
      if projection is None or len(projection) <= 0 :
         projection = "SIN" 
      postfix = ( "_%.2fMHz_BeamCorr%s.txt" % (options.frequency_mhz,projection.upper()) ) 
      outfile=options.time_azh_file.replace(".txt",postfix)
      beam_correct_flux( dt, tm, az, el, ra, dec, flux, cnt, frequency_mhz=options.frequency_mhz, outfile=outfile, pol=options.polarisation, projection=options.projection )

   elif options.lightcurve_file is not None :
      print("Reading lightcurve file %s ..." % (options.lightcurve_file))
      
      dt_index  = 0
      tm_index  = 0 
      az_index  = 9
      el_index  = 10
      ra_index  = -1
      dec_index = -1
      flux_index = 1 
      
      if options.dt_index is not None :
         dt_index = options.dt_index

      if options.tm_index is not None :
         dt_index = options.tm_index

      if options.az_index is not None :
         az_index = options.az_index
      
      if options.el_index is not None :
         el_index = options.el_index

      if options.ra_index is not None :
         ra_index = options.ra_index

      if options.dec_index is not None :
         dec_index = options.dec_index

      if options.flux_index is not None :
         flux_index = options.flux_index
      
      (dt,tm,az,el,ra,dec,flux,cnt) = read_time_azh_file( options.lightcurve_file , dt_index=dt_index, tm_index=tm_index,   az_index=az_index, el_index=el_index, 
                                                                                    ra_index=ra_index, dec_index=dec_index, flux_index=flux_index,
                                                                                    default_ra=options.ra_deg, default_dec=options.dec_deg )
      projection = options.projection
      if projection is None or len(projection) <= 0 :
         projection = "SIN" 
      postfix = ( "_%.2fMHz_BeamCorr%s.txt" % (options.frequency_mhz,projection.upper()) ) 
      outfile=options.lightcurve_file.replace(".txt",postfix)
      beam_correct_flux( dt, tm, az, el, ra, dec, flux, cnt, frequency_mhz=options.frequency_mhz, outfile=outfile, pol=options.polarisation, projection=options.projection )
      
      
   elif options.fits2beam is not None :
      print("DEBUG : fits2beam for file %s" % (options.fits2beam))      
      fits2beam( options.fits2beam, options )
      
   elif options.matlab2fits is not None :   
      matlab2sin( options.matlab2fits, x_size=512 )
      
   elif options.infile_hdf5list is not None : 
      save_beam_on_sun_file( options )
      
   else :
      az = options.azim_deg
      za = options.za_deg
      freq_mhz = options.frequency_mhz
      
      if options.beam_on_sun :
         if options.unix_time is None or options.unix_time <= 0 :
            print("ERROR : calculation of beam on Sun is required, but time is not specified -> cannot continue, please provide --uxtime or --unix_time parameter value")
            os.sys.exit(-1)
      
         print("INFO : calculation of beam on Sun is required")         
         ( ra, dec, az, alt, za ) = sun_position( options.unix_time )
         print("INFO : sun position calculated (RA,DEC) = (%.4f,%.4f) [deg] , (AZ,ELEV,ZA) = (%.4f,%.4f,%.4f) [deg]" % (ra, dec, az, alt, za))
      
      beam_x = 0
      try : 
         beam_x = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='X', projection=options.projection, station_name=options.station_name )   
      except : 
         print("ERROR : beam-X cannot be calculated")

      beam_y = 0         
      try :          
         beam_y = get_fits_beam( numpy.array([[az]]) , numpy.array([[za]]) , freq_mhz, polarisation='Y', projection=options.projection, station_name=options.station_name )   
      except :
         print("ERROR : beam-Y cannot be calculated")
         
      print("BEAM_X = %.4f , BEAM_Y = %.4f " % (beam_x,beam_y))
      
      beam_x_2 = [[0]]
      try : 
         beam_x_2 = get_fits_beam_multi( numpy.array([[az*(math.pi/180.00)]]) , numpy.array([[za*(math.pi/180.00)]]) , freq_mhz, polarisation='X' , station_name=options.station_name, projection=options.projection )
      except :
         print("ERROR : beam-X cannot be calculated using get_fits_beam_multi")
         
      beam_y_2 = [[0]] 
      try :    
         beam_y_2 = get_fits_beam_multi( numpy.array([[az*(math.pi/180.00)]]) , numpy.array([[za*(math.pi/180.00)]]) , freq_mhz, polarisation='Y' , station_name=options.station_name, projection=options.projection )
      except :
         print("ERROR : beam-Y cannot be calculated using get_fits_beam_multi")
         
      print("DEBUG2: BEAM_X = %.4f , BEAM_Y = %.4f ( SQUARED = %.4f , %.4f )" % (beam_x_2[0][0],beam_y_2[0][0],(beam_x_2[0][0])**2,(beam_y_2[0][0])**2))

      
      
      