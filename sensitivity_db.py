from __future__ import print_function
# SCRIPT for getting station sensitivity from the sensitivity database (SQLITE3 or PostgreSQL) :
# example commands :
#  1/ A/T vs. frequency at a particular pointing direction :
# 
#    python ./sensitivity_db.py --azim_deg=0 --za_deg=0 --lst=15.4 --out_file="Azim0_Za0_Lst15.4hours" --do_plot 
#    python ./sensitivity_db.py --azim_deg=0 --za_deg=30 --lst=15.4 --out_file="Azim0_Za0_Lst15.4hours" --do_plot
#    Verify results against files in templates/ Azim0_Za0_Lst15.4hours_X.txt and Azim0_Za0_Lst15.4hours_Y.txt
#       cat nimbus5/za30-60deg/1262476817_az0_za30_eda_sensitivity_YY.txt
#       cat comparison_AAVS2_EDA2/EDA2/Galaxy_transit/1262487044_az0_za30_eda_sensitivity_YY.txt
#       cat comparison_AAVS2_EDA2/EDA2/Galaxy_nadir/1262530097_az0_za30_eda_sensitivity_XX.txt
#
#     pointing specified with (ra,dec) , example for 3C444 :
#        python ./sensitivity_db.py --ra=333.60724995 --dec=-17.02661111 --lst=19.33007231 --out_file=3c444_test --do_plot
# 
# 
#  2/ Create map of sensitivity for lst=15.4 hours and freq = 154.88 MHz :
#     python ./sensitivity_db.py --freq_mhz=154.88 --lst_hours=15.4
#     python ./sensitivity_db.py --freq_mhz=154.88 --lst_hours=15.4 --save_text
#     python ./sensitivity_db.py --freq_mhz=154.88 --lst_hours=2.00 --do_plot
#     Galactic transit : python ./sensitivity_db.py --freq_mhz=154.88 --lst_hours=17.76 --do_plot
#     Galsctic nadir   : python ./sensitivity_db.py --freq_mhz=154.88 --lst_hours=5.73 --do_plot
# 
#  3/ A/T vs. time :
#     a/ A/T vs. unix time / UTC at a particular pointing direction and frequency at zenith 
#        python ./sensitivity_db.py --freq_mhz=154.88 --unixtime_start=1582597729  --interval=86400 --azim_deg=0 --za_deg=0 --do_plot
# 
#        pointing specified with (ra,dec) , example for 3C444 :
#        python ./sensitivity_db.py --freq_mhz=160.00 --unixtime_start=1636777350 --interval=86400 --ra=333.60724995 --dec=-17.02661111 --do_plot
#        
#
#     b/ A/T vs. lst at a particular pointing direction and frequency at zenith
#        - for pointing specified by (AZ,EL) :          
#          python ./sensitivity_db.py --freq_mhz=154.88 --lst_start=0.00 --lst_end=24.00 --azim_deg=0 --za_deg=0 --do_plot      
#  
#        - for pointing specified by (RA,DEC) in degrees, example for 3C444 :
#          python ./sensitivity_db.py --freq_mhz=160.00 --lst_start=0.00 --lst_end=24.00 --ra=333.60724995 --dec=-17.02661111 --do_plot
# 
#    Sensitivity map :
#   https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.SmoothSphereBivariateSpline.html
# import sensitivity_db
# (azim_x,za_x,aot_x,sefd_x,azim_y,za_y,aot_y,sefd_y) = sensitivity_db.get_sensitivity_map( 160.0000, 15.4 )
# sensitivity_db.plot_sensitivity_map( azim_x,za_x,aot_x )
# or try :
# from scipy.interpolate import SmoothSphereBivariateSpline 
# import numpy
# import numpy as np
# import matplotlib.pyplot as plt
# azim_x=azim_x*(numpy.pi/180.0)
# za_x=za_x*(numpy.pi/180.0)
# lut = SmoothSphereBivariateSpline( za_x, azim_x,aot_x , s=3.5)
# fine_lats = np.linspace(0., np.pi, 90)
# fine_lons = np.linspace(0., 2 * np.pi, 360)
# data_smth = lut(fine_lats, fine_lons)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1 = fig.add_subplot(111,polar=True)
# ax1.imshow(data_smth, interpolation='nearest')
# plt.show()
# or : Google : "matplotlib polar heatmap" : https://stackoverflow.com/questions/36513312/polar-heatmaps-in-python


import sqlite3
from optparse import OptionParser

# plotting :
from pylab import *
import numpy
import math
import matplotlib.pyplot as plt
import matplotlib.dates as md  # to convert unix time to date on X-axis 
import os
from scipy.interpolate import interp1d # for interpolation 

import matplotlib
matplotlib.rc('xtick', labelsize=25) 
matplotlib.rc('ytick', labelsize=25) 

# local packages :
import beam_tools
import fits_beam

# time operations :
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from datetime import datetime
MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)
# unixtime -> string -> gps :
#   t_ux=datetime.utcfromtimestamp(1582274015).strftime('%Y-%m-%d %H:%M:%S')
#   t_now = Time( t_ux, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))
#   t_now.gps
#   t_now.unix
# 

# code version
code_version = 1.00
code_version_line = "# Data generated with code version 1.00\n"

debug_level = 2

# plot SKA-Low requirements 
plot_requirements = True
import lfaa_requirements

plot_braun2019 = True
import sensitivity_braun2019

# path where images are saved, by default local directory :
# save_output_path = "./"

# script for quering SQLITE3 or PostgreSQL databases for sensitivity values :
# HELP : python : https://www.sqlitetutorial.net/sqlite-python/ , https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/
#        Plot vs. time : /home/msok/ska/aavs/aavs0.5/trunk/analysis/MWAdas/scripts$ plot_delay.py
# 
# Google : "python interpolation on sphere"
#        https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.interpolate.SmoothSphereBivariateSpline.html
#        https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.LSQSphereBivariateSpline.html
#        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.SmoothSphereBivariateSpline.html


web_interface_initialised = False
try:
##   from io import StringIO
   from io import BytesIO
except ImportError:
   print("WARNING : import ( from io import StringIO ) failed")

web_interface_initialised = False 
def init_web_interface() :
   global web_interface_initialised

   if not web_interface_initialised :
#      import cStringIO 
#      from io import StringIO
#      import base64 # Use this in Python 3x
      web_interface_initialised = True

def create_connection_sqlite3( db_file ):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        print("Connecting to database %s" % (db_file))
        conn = sqlite3.connect(db_file)
    except sqlite3.Error as e:
        print(e)        
        print("ERROR : could not connect to database %s" % (db_file))

 
    return conn

def calc_anglular_distance_degrees( azim1_deg, za1_deg, azim2_deg , za2_deg ) :
    azim1_rad = azim1_deg * (math.pi/180.00)
    elev1_rad = ( 90.00 - za1_deg ) * (math.pi/180.00)

    azim2_rad = azim2_deg * (math.pi/180.00)
    elev2_rad = ( 90.00 - za2_deg ) * (math.pi/180.00)

    cos_value = math.sin(elev1_rad)*math.sin(elev2_rad) + math.cos(elev1_rad)*math.cos(elev2_rad)*math.cos( azim1_rad - azim2_rad )
    dist_value_deg = math.acos( cos_value )*(180.00/math.pi)
    
    return dist_value_deg

def aotxy2aoti( aot_x , aot_y ):
   k_b = 1380.00
   sefd_x = (2*k_b)/aot_x
   sefd_y = (2*k_b)/aot_y
   sefd_i = 0.5*math.sqrt( sefd_x*sefd_x + sefd_y*sefd_y )
   aot_i = (2*k_b)/sefd_i
   
   return (aot_i,sefd_i)


def ha2range_deg( ha ) :
  if ha < 0.00 :
     ha = ha + 360.00;
  elif ha > 360.00 :
     ha = ha - 360.00 
     
  return ha 

# based on /opt/pi/dev/pisys/daq/src/ccd/ccdastro/AstroCCD.cpp void AstroCCD::calculateHourAngleBase( double ra, double dec, double& sid_time, double& ha )
# 
#  calculates (AZIM,ELEVATION) from (RA,DEC) using LST 
# 
def radec2azim( ra_deg, dec_deg, lst_hours, geo_lat=-26.70331944444445, debug=True, astro_azim=True ) : # default GEO_LAT is MRO 
   DEG2RAD = (math.pi/180.00)
   RAD2DEG = (180.00/math.pi)

   ha_deg = lst_hours*15.00 - ra_deg
   ha_deg = ha2range_deg( ha_deg )
   if debug : 
      print("DEBUG : ha := %.4f [deg]" % (ha_deg))
         
   sin_alt = math.sin( geo_lat*DEG2RAD ) * math.sin( dec_deg*DEG2RAD )  + math.cos( geo_lat*DEG2RAD ) * math.cos( dec_deg*DEG2RAD ) * math.cos( ha_deg*DEG2RAD )
   alt_rad = math.asin( sin_alt )
   alt_deg = alt_rad * RAD2DEG
   
   up = ( math.cos( dec_deg*DEG2RAD ) * math.sin( ha_deg*DEG2RAD ) ) 
   bottom =  ( math.sin( geo_lat*DEG2RAD )*math.cos( dec_deg*DEG2RAD )*math.cos( ha_deg*DEG2RAD ) - math.cos( geo_lat*DEG2RAD )*math.sin( dec_deg*DEG2RAD ) ) 
   azim_rad = math.atan2( up, bottom )
   azim_deg = azim_rad*RAD2DEG
   
   if astro_azim : 
      azim_deg = (180.00 + azim_deg)
   
   print("(Azim,alt) = (%.4f,%.4f) [deg]" % (azim_deg,alt_deg))
   
   return ( azim_deg, alt_deg, ha_deg, ra_deg, dec_deg, lst_hours, geo_lat )


def read_text_file( filename , do_fit=True ) : 
   print("read_data(%s) ..." % (filename))

   if not os.path.exists( filename ) :
      print("ERROR : could not read satellite info file %s" % (filename))
      return (None,None,None)
   
   file=open(filename,'r')

   # reads the entire file into a list of strings variable data :
   data=file.readlines()
   # print data
   
   freq_mhz = []
   t_rcv    = []

   for line in data : 
      line=line.rstrip()
      words = line.split(' ')

      if line[0] == '#' or line[0]=='\n' or len(line) <= 0 or len(words)<2 :
         continue

#      print "line = %s , words = %d" % (line,len(words))

      if line[0] != "#" :
         f  = float( words[0+0] )
         t = float( words[1+0] )
         
         freq_mhz.append( f )
         t_rcv.append( t )


   freq_mhz = numpy.array( freq_mhz )
   t_rcv    = numpy.array( t_rcv )

   fit_func = None
   if do_fit :
      fit_func = interp1d( freq_mhz, t_rcv, kind='cubic')  
      #  ret_val = lna_gain_db_cubic( freq_mhz )
         
         
   return (freq_mhz,t_rcv,fit_func)

# calculate Stokes I using an APPROXIMATE formula, but this is best we can do see Sutijno, Sokolowski et al (2020) for the improved one :
def calc_sefd_i( out_sefd_x, out_freq_x, out_sefd_y, out_freq_y, out_file=None ) :
    out_freq_i = None
    out_sefd_i = []
    out_aot_i = []
    
    if out_freq_x is not None :
       out_freq_i = []

    # calculate SEFD_I - if possible (same sizes of arrays):
    if len(out_sefd_x) == len(out_sefd_y) :
       for i in range(0,len(out_sefd_x)) :
          sefd_x = out_sefd_x[i]
          sefd_y = out_sefd_y[i]
          
          sefd_i = 0.5*math.sqrt( sefd_x*sefd_x + sefd_y*sefd_y )

          if out_freq_x is not None :
             out_freq_i.append( out_freq_x[i] )
          out_sefd_i.append( sefd_i )

          aot_i = 0
          if not numpy.isnan(sefd_i) :
             aot_i = (2*1380.00)/sefd_i
             
          out_aot_i.append( aot_i )          
          
          if out_file is not None :
             line = "%.4f %.4f %.8f\n" % (out_freq_x[i],out_freq_y[i],aot_i)
             out_file.write( line )
          
#          print("DEBUG : %.4f %.4f -> %.4f -> %.4f" % (sefd_x,sefd_y,sefd_i,aot_i))
    else :
       print("ERROR : len(out_aot_x) != len(out_aot_y) ( %d != %d )" % (len(out_aot_x),len(out_aot_y)))


    return (out_freq_i,out_sefd_i,out_aot_i)

# get sensitivity vs. frequency for given pointing direction [degrees] and lst [hours]:
def get_sensitivity_azzalst( az_deg , za_deg , lst_hours , 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00,
                             receiver_temp_file=None, 
                             ra_deg=None, dec_deg=None ) :
        
    global debug_level
    
    print("DEBUG : get_sensitivity_azzalst")
    
    if ra_deg is not None and dec_deg is not None :
       # def radec2azim( ra_deg, dec_deg, lst_hours, geo_lat=-26.70331944444445, debug=True, astro_azim=True ) : # default GEO_LAT is MRO 
       # return ( azim_deg, alt_deg, ha_deg, ra_deg, dec_deg, lst_hours, geo_lat )
       
       (az_deg, alt_deg, ha_deg, ra_deg2, dec_deg2, lst_hours, geo_lat ) = radec2azim( ra_deg, dec_deg, lst_hours, geo_lat=MWA_POS.lat.value )
       za_deg = (90.00 - alt_deg )
       print("DEBUG : calculated horizontal coordinates (az,za,ha) = (%.4f,%.4f,%.4f) [deg] from (ra,dec) = (%.4f,%.4f) [deg] at LST = %.4f [h]" % (az_deg, za_deg, ha_deg, ra_deg, dec_deg, lst_hours)) 
    
    # 
    file_freq_mhz = None
    file_trcv     = None
    file_trcv_func = None
    if receiver_temp_file is not None :
       (file_freq_mhz,file_trcv,file_trcv_func) = read_text_file( receiver_temp_file )
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
    
    if conn is None :
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None)

  
    # select closest LST :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT MIN(ABS(lst-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f or %.4f<0.1)" %  (lst_hours,lst_hours,db_lst_resolution, za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg )
    print("DEBUG SQL1 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_lst_distance = None
    for row in rows:
       if row[0] is not None :
          min_lst_distance = float( row[0] )
       else :
          print("ERROR : not data close to lst = %.2f h at (az,za) = (%.4f,%.4f) [deg] in the database" % (lst_hours,az_deg,za_deg))

    if min_lst_distance is None :
       print("ERROR no records in the database exist closer than %.4f hours in LST at (az,za) = (%.4f,%.4f) [deg]" % (db_lst_resolution,az_deg,za_deg))
       return (None,None,None,None,None,None)

    print("Best LST in database is closer than %.8f [hours]" % (min_lst_distance))
    min_lst_distance = min_lst_distance + 0.01
    
    if min_lst_distance is None :
       print("ERROR no records in the database exist closer than %.4f hours in LST in the direction of (azim,za) = (%.4f,%.4f) [deg]" % (db_lst_resolution,az_deg,za_deg))
       return (None,None,None,None,None,None)
 
    # get requested data :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f or %.4f<0.1)" %  (lst_hours,(min_lst_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg , za_deg )
    print("DEBUG SQL2 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    
    # find closest pointing direction (WARNING : no cos and sin in SQLITE3 ):
    min_angular_distance_deg = 1e6
    closest_gridpoint_za_deg = -10000.00
    closest_gridpoint_az_deg = -10000.00
    for row in rows:
       azim_deg_db = float( row[1] )
       za_deg_db   = float( row[2] )
       
       # CalcDistRADEC( ra1, dec1, ra2, dec2 )
       dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
       
       if dist_value_deg < min_angular_distance_deg :
          min_angular_distance_deg = dist_value_deg
          closest_gridpoint_za_deg = za_deg_db
          closest_gridpoint_az_deg = azim_deg_db

    if min_angular_distance_deg > 1000 :
       print("ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg))
       return (None,None,None,None,None,None)
       
       
    print("Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg))
     
    out_freq_x = []
    out_aot_x  = []
    out_sefd_x = []
    
    out_freq_y = []
    out_aot_y  = []
    out_sefd_y = []

    for row in rows:
        if debug_level > 0 : 
           print(row)
        
        id = int( row[0] )
        azim_deg_db = float( row[1] )
        za_deg_db   = float( row[2] )
        freq_mhz    = float( row[3] )
        pol         = row[4]
        lst_db      = float( row[5] )
        unixtime    = float( row[6] )
        gpstime     = float( row[7] )
        aot         = float( row[8] )
        sefd        = numpy.NaN
        if aot > 0 :
           sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        if file_trcv_func is not None :
           # in case trcv is provided in the external text file (like a config file) :
           t_rcv = file_trcv_func( freq_mhz )
           aot = a_eff / ( t_ant + t_rcv )                
           sefd = 0.00
           if aot > 0 :
              sefd        = (2*1380.00)/aot
        
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db, az_deg, za_deg )
              
        print("TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,ang_distance_deg,az_deg,za_deg))
        
        
        
        if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
            if pol == "X" :
               out_freq_x.append( freq_mhz )
               out_aot_x.append( aot )
               out_sefd_x.append( sefd ) 
            elif pol == "Y" :
               out_freq_y.append( freq_mhz )
               out_aot_y.append( aot )
               out_sefd_y.append( sefd )            
            else :
               print("ERROR : unknown polarisation = %s" % (pol))
        else :
           print("\t\tDB Record ignored due to angular distance %.4f [deg] > limit = %.4f [deg]" % (ang_distance_deg,(min_angular_distance_deg+0.01)))
 

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_freq_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_freq_x, out_sefd_y, out_freq_y )

    return ( numpy.array(out_freq_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_freq_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_freq_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )


def unixtime2lst( unixtime ) :
    utc_str = datetime.utcfromtimestamp( unixtime ).strftime('%Y-%m-%d %H:%M:%S')
    t_utc = Time( utc_str, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))

    lst_hours = t_utc.sidereal_time('apparent').value
    
    print("DEBUG : unixtime2lst( %d ) = %.4f [hours]" % (unixtime,lst_hours))
    
    return lst_hours


# get sensitivity vs. unix time for a specific polarisation :
def get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=5.00, 
                             receiver_temperature=None,
                             ra_deg=None, dec_deg=None ) :

    print("DEBUG : get_sensitivity_timerange_single_pol")

    horizontal_coord = True
    if ra_deg is not None and dec_deg is not None :
       horizontal_coord = False
       
    out_uxtime = []
    out_aot    = []
    out_sefd   = []
    
    out_f = None
    if pol == "X" or pol == "Y" :
       tmp_filename = "test_%s.txt" % (pol)
       out_f = open( tmp_filename ,"w")
    

    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
    
    if conn is None :
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None, None, None)

  
    # get requested data :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f or %.4f<0.1) AND polarisation='%s'" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol )
    print("DEBUG SQL2 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    
# TODO : needs to be done per polarisation to easily identify best_row_id and use it - otherwise to complicted. Create a function     
    
    # find closest pointing direction at any LST (WARNING : no cos and sin in SQLITE3 ):
    min_angular_distance_deg = db_ang_res_deg

# MS : 2021-11-29 : calculating min_angular_distance_deg for each timestamps to make it working with A/T vs. time for RADEC     
#    if ra_deg is not None and dec_deg is not None :
#       print("INFO : pointing specified by (ra,dec) = (%.4f,%.4f) [deg] -> min_angular_distance_deg = %.4f [deg]" % (ra_deg,dec_deg,min_angular_distance_deg))
#    else :
#       print("INFO : pointing specified by (azim,za) = (%.4f,%.4f) [deg] -> looking for min_angular_distance_deg ..." % (az_deg,za_deg))
#       closest_gridpoint_za_deg = -10000.00
#       closest_gridpoint_az_deg = -10000.00
#       for row in rows:
#          azim_deg_db = float( row[1] )
#          za_deg_db   = float( row[2] )
#       
#          # CalcDistRADEC( ra1, dec1, ra2, dec2 )
#          dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
#       
#          if dist_value_deg < min_angular_distance_deg :
#             min_angular_distance_deg = dist_value_deg
#             closest_gridpoint_za_deg = za_deg_db
#             closest_gridpoint_az_deg = azim_deg_db
#             
#          print("Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg))   
#          
#       print("INFO : pointing specified by (azim,elev) = (%.4f,%.4f) [deg] -> min_angular_distance_deg = %.4f [deg]" % (az_deg,za_deg,min_angular_distance_deg))   
#
#    if min_angular_distance_deg > 1000 :
#       print("ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg))
#       return (None,None,None,None,None,None)
       
       

# MS : 2021-11-29 - modified moved to inside the loop over timetamps :
# when used_db_record_id is defined before the loop over timestamps there will not be flat A/T parts ------ when same A/T record is used for multiple timestamps 
    used_db_record_id = {}
    for unixtime in range( int(ux_start), int(ux_end)+1, time_step ) :
       # for a given timestamp a record used only once :
#       used_db_record_id = []
 
       lst_hours = unixtime2lst( unixtime )
       
       # if (RA,DEC) are specified -> use them to calculate AZIM,ELEV for the current LST :
       if ra_deg is not None and dec_deg is not None :
          (az_deg, alt_deg, ha_deg, ra_deg2, dec_deg2, lst_hours, geo_lat ) = radec2azim( ra_deg, dec_deg, lst_hours, geo_lat=MWA_POS.lat.value )
          za_deg = (90.00 - alt_deg )
          print("DEBUG : calculated horizontal coordinates (az,za,ha) = (%.4f,%.4f,%.4f) [deg] from (ra,dec) = (%.4f,%.4f) [deg] at LST = %.4f [h]" % (az_deg, za_deg, ha_deg, ra_deg, dec_deg, lst_hours)) 

       
#       print("DEBUG : ux=%.2f -> lst=%.2f" % (unixtime,lst_hours))
    
       cur = conn.cursor()       
       # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
       szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,a_eff,t_sys,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE lst>=%.4f AND lst<=%.4f AND frequency_mhz>=%.4f AND frequency_mhz<=%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND polarisation='%s' ORDER BY za_deg,azim_deg,lst" %  \
                 ( (lst_hours-db_lst_resolution),(lst_hours+db_lst_resolution),\
                   (freq_mhz-(db_freq_resolution_mhz+0.01)), (freq_mhz+(db_freq_resolution_mhz+0.01)),\
                   za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol\
                  )
       if debug_level >= 2 :
          print("DEBUG SQL get_sensitivity_timerange_single_pol : %s" % (szSQL))
       cur.execute( szSQL )
       rows = cur.fetchall()
       
       min_lst_distance = db_lst_resolution
       best_rec_id = -1
       best_rec_id_list = []
       if len(rows) > 1 :               
           # TODO : should only be one row for both polarisations (select the closest in time)
           for row in rows :
              id = int( row[0] )
              azim_deg_db = float( row[1] )
              za_deg_db   = float( row[2] )
              freq_mhz    = float( row[3] )
              pol         = row[4]
              lst_db      = float( row[5] )

              if math.fabs( lst_db - lst_hours ) < min_lst_distance :
                 min_lst_distance = math.fabs( lst_db - lst_hours )
                 best_rec_id = id
                 
                 if debug_level >= 2 :
                    print("DEBUG : lst_db = %.4f h, lst requested = %.4f h -> min_lst_distance = %.4f h (azim,za) = (%.4f,%.4f) [deg] -> best_rec_id = %d" % (lst_db,lst_hours,min_lst_distance,azim_deg_db,za_deg_db,best_rec_id))
              
              dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
              if dist_value_deg < min_angular_distance_deg :
                 min_angular_distance_deg = dist_value_deg

           # Repeat the loop again because there may be records with the same exact LST which we should use because they have the smallest angular distance from out pointing direction:
           for row in rows :
              id = int( row[0] )
              azim_deg_db = float( row[1] )
              za_deg_db   = float( row[2] )
              freq_mhz    = float( row[3] )
              pol         = row[4]
              lst_db      = float( row[5] )

              if math.fabs( lst_db - lst_hours ) <= (min_lst_distance+0.1) :
                 best_rec_id_list.append( id )

           
                 
       if debug_level >= 2 or min_lst_distance < db_lst_resolution :
          print("DEBUG : closest in time is min_lst_diff = %.4f [h] (id = %d)" % (min_lst_distance,best_rec_id))


       min_angular_distance_deg = db_ang_res_deg
       best_uxtime = -1
       best_aot = -1
       best_sefd = -1
       best_id = -1
       best_lst_db = -1

       last_lst = -1000        
       for row in rows :
           if debug_level > 0 : 
              print(row)
        
           id = int( row[0] )
           lst_db      = float( row[5] )

           # if duplicates not skipped -> introduces "plataus" of the same value for several times :           
           if id in used_db_record_id.keys() :
              # skip already used records :
              continue 
#              if horizontal_coord or fabs( lst_db - lst_hours) >= used_db_record_id[id] :
#                 continue
#              else :
#                 if debug_level >= 2 :
#                    print("DEBUG : record id = %d (lst = %.4f) already used, but allowing because the timestamp is closer to %.4f" % (id,lst_db,lst_hours)) 
           
           azim_deg_db = float( row[1] )
           za_deg_db   = float( row[2] )
           freq_mhz    = float( row[3] )
           pol         = row[4]
           lst_db      = float( row[5] )
           unixtime_db = float( row[6] ) # when filled 
           gpstime     = float( row[7] )
           aot         = float( row[8] )
           sefd        = numpy.NaN
           if aot > 0 :
              sefd        = (2*1380.00)/aot
           a_eff       = float( row[9] )
           t_sys       = float( row[10] )
           t_rcv       = float( row[11] )
           t_ant       = float( row[12] )
           array_type  = int( row[13] )
           timestamp   = row[14]
           creator     = row[15]
           code_version = row[16]
           
           dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
           
           if receiver_temperature is not None and receiver_temperature >= 0.00 :
              # in case trcv is provided in the external text file (like a config file) :
              t_rcv = receiver_temperature
              aot = a_eff / ( t_ant + t_rcv )
              sefd        = (2*1380.00)/aot
           
           print("TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg] (id = %d vs. best_rec_id = %d)" % (id,freq_mhz,azim_deg_db,za_deg_db,dist_value_deg,az_deg,za_deg,id,best_rec_id))
                
           if best_rec_id < 0 or id in best_rec_id_list : 
              ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg ) 
        
              # if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
              if ang_distance_deg <= min_angular_distance_deg :
                  if debug_level >= 4 :
                     is_in_list = False
                     if id in used_db_record_id.keys() :
                        is_in_list = True
                     print("DEBUG_FINAL id = %d : %.4f %.4f %.4f , list count = %d , is_in_list = %s" % (id,unixtime,aot,sefd,len(used_db_record_id),is_in_list))
                     
#                  out_uxtime.append( unixtime )
#                  out_aot.append( aot )
#                  out_sefd.append( sefd )                   
#                  used_db_record_id.append( id ) 
                  best_uxtime = unixtime
                  best_aot = aot
                  best_sefd = sefd
                  best_id = id   
                  best_lst_db = lst_db
                  if debug_level >= 2 : 
                     print("\t\tDEBUG : updated best record ID = %d" % (best_id))
              else :
                 print("\t\tDB Record ignored due to angular distance = %.4f [deg] larger than limit of %.4f [deg]" % (ang_distance_deg,(min_angular_distance_deg+0.01)))
           else :
              print("\t\tRecord %d not checked due to best_rec_id = %d >= 0 or ID not in list !" % (id,best_rec_id))
             
           if out_f is not None :
              line = "%.4f %.4f\n" % (unixtime,t_ant) 
              out_f.write( line )
           
           last_lst = lst_db
 
       if best_id > 0 :
          print("DEBUG pol=%s : UXTIME = %.2f closest in time and space is ID = %d , aot = %.4f [m^2/K]" % (pol,best_uxtime,best_id,best_aot))
          out_uxtime.append( best_uxtime )
          out_aot.append( best_aot )
          out_sefd.append( best_sefd )
          used_db_record_id[ best_id ] = 1 # math.fabs( best_lst_db - lst_hours ) #
 
    if out_f is not None :
       out_f.close()
              
    if debug_level >= 2 :
       print("DEBUG information on sensitivity :")
       for i in range(0,len(out_uxtime)) :
          print("%.2f %.4f %.4f" % (out_uxtime[i],out_aot[i],out_sefd[i]))   
          
    print("DEBUG : get_sensitivity_timerange_single_pol pol = %s returns %d / %d / %d records" % (pol,len(out_uxtime),len(out_aot),len(out_sefd)))         

    return (out_uxtime, out_aot, out_sefd)

# get sensitivity vs. LST for a specific polarisation :
def get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None ) :

    print("DEBUG : get_sensitivity_lstrange_single_pol")

    out_lst = []
    out_aot    = []
    out_sefd   = []

    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
    
    if conn is None :
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None)              

    # select closest FREQUENCT :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT MIN(ABS(frequency_mhz-%.4f)) FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol, lst_start, lst_end )
    print("DEBUG SQL1 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_freq_distance = None
    for row in rows:
       if row[0] is not None :
          min_freq_distance = float( row[0] )
       else :
          print("ERROR : not data close to frequency = %.2f [MHz] at (az,za) = (%.4f,%.4f) [deg] in the database" % (freq_mhz,az_deg,za_deg))

    print("DEBUG : minimum frequency distance = %.2f MHz" % (min_freq_distance))
  
    # get requested data :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<=%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,(min_freq_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol, lst_start, lst_end )
    print("DEBUG SQL2 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    
# TODO : needs to be done per polarisation to easily identify best_row_id and use it - otherwise to complicted. Create a function     
    
    # find closest pointing direction at any LST (WARNING : no cos and sin in SQLITE3 ):
    min_angular_distance_deg = db_ang_res_deg
    closest_gridpoint_za_deg = -10000.00
    closest_gridpoint_az_deg = -10000.00
    for row in rows:
       azim_deg_db = float( row[1] )
       za_deg_db   = float( row[2] )
       
       # CalcDistRADEC( ra1, dec1, ra2, dec2 )
       dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
       
       if dist_value_deg < min_angular_distance_deg :
          min_angular_distance_deg = dist_value_deg
          closest_gridpoint_za_deg = za_deg_db
          closest_gridpoint_az_deg = azim_deg_db

    if min_angular_distance_deg > 1000 :
       print("ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg))
       return (None,None,None,None,None,None)
       
       
    print("Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg))


    
    cur = conn.cursor()       
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND polarisation='%s' AND lst>%.2f AND lst<%.2f ORDER BY LST ASC" %  (freq_mhz,(min_freq_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol, lst_start, lst_end )
    if debug_level >= 2 :
       print("DEBUG SQL get_sensitivity_timerange_single_pol : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
       
    for row in rows :
        if debug_level > 0 : 
           print(row)
        
        id = int( row[0] )
        azim_deg_db = float( row[1] )
        za_deg_db   = float( row[2] )
        freq_mhz    = float( row[3] )
        pol         = row[4]
        lst_db      = float( row[5] )
        unixtime    = float( row[6] )
        gpstime     = float( row[7] )
        aot         = float( row[8] )
        if aot != 0 :
           sefd        = (2*1380.00)/aot
        else :
           sefd     = 1e20
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        if receiver_temperature is not None and receiver_temperature >= 0.00 :
           t_rcv = receiver_temperature
           aot = a_eff / ( t_ant + t_rcv )
           sefd        = (2*1380.00)/aot
           
        print("TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,min_angular_distance_deg,az_deg,za_deg))
                
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg ) 
        
        if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
            out_lst.append( lst_db )
            out_aot.append( aot )
            out_sefd.append( sefd ) 
        else :
            print("\t\tDB Record (azim,za) = (%.4f,%.4f) [deg] ignored due to angular distance too large ( %.4f deg )" % (azim_deg_db,za_deg_db,ang_distance_deg))
 
    return (out_lst, out_aot, out_sefd)


# get sensitivity vs. LST for a specific polarisation :
def get_sensitivity_radec_lstrange_single_pol( ra_deg , dec_deg , freq_mhz, lst_start, lst_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None ) :

    print("DEBUG : get_sensitivity_lstrange_single_pol")

    out_lst = []
    out_aot    = []
    out_sefd   = []

    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
    
    if conn is None :
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None)              

    # select closest FREQUENCY :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT MIN(ABS(frequency_mhz-%.4f)) FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,freq_mhz,(db_freq_resolution_mhz+0.01), pol, lst_start, lst_end )
    print("DEBUG SQL1 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_freq_distance = None
    for row in rows:
       if row[0] is not None :
          min_freq_distance = float( row[0] )
       else :
          print("ERROR : not data close to frequency = %.2f [MHz] in LST range %.4f - %.4f [h] the database" % (freq_mhz,lst_start, lst_end))

    print("DEBUG : minimum frequency distance = %.2f MHz" % (min_freq_distance))
  
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    different_lsts = []
    szSQL = "SELECT DISTINCT(lst) FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<=%.4f AND polarisation='%s' AND lst>=%.2f AND lst<=%.2f" %  (freq_mhz,(min_freq_distance+0.01), pol, lst_start, lst_end )
    print("DEBUG SQL2 : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()

    print("Different LSTs:")    
    for row in rows:
       lst_value = float( row[0] )
       different_lsts.append( lst_value )
       print("\t%.4f" % (lst_value))

    # loop over LSTs in the range : 
    for lst in different_lsts :         
       cur = conn.cursor()       
       
       (az_deg, alt_deg, ha_deg, ra_deg2, dec_deg2, lst_hours, geo_lat ) = radec2azim( ra_deg, dec_deg, lst, geo_lat=MWA_POS.lat.value )
       za_deg = (90.00 - alt_deg )
       
       # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
       szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND polarisation='%s' AND ABS(lst-%.8f)<0.0001 ORDER BY ((azim_deg-%.8f)*(azim_deg-%.8f)+(za_deg-%.8f)*(za_deg-%.8f)) ASC" %  (freq_mhz,(min_freq_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, pol, lst, az_deg, az_deg, za_deg, za_deg )
       if debug_level >= 2 :
          print("DEBUG SQL get_sensitivity_radec_lstrange_single_pol : %s" % (szSQL))
       cur.execute( szSQL )
       rows = cur.fetchall()
       
       min_angular_distance_deg = db_ang_res_deg
       
       best_id = -1
       best_lst = -1
       best_aot = -1
       best_sefd = -1
       
       for row in rows :
           if debug_level > 0 : 
              print(row)
        
           id = int( row[0] )
           azim_deg_db = float( row[1] )
           za_deg_db   = float( row[2] )
           freq_mhz    = float( row[3] )
           pol         = row[4]
           lst_db      = float( row[5] )
           unixtime    = float( row[6] )
           gpstime     = float( row[7] )
           aot         = float( row[8] )
           if aot != 0 :
              sefd        = (2*1380.00)/aot
           else :
              sefd     = 1e20
           a_eff       = float( row[9] )
           t_rcv       = float( row[10] )
           t_ant       = float( row[11] )
           array_type  = int( row[12] )
           timestamp   = row[13]
           creator     = row[14]
           code_version = row[15]
        
           if receiver_temperature is not None and receiver_temperature >= 0.00 :
              t_rcv = receiver_temperature
              aot = a_eff / ( t_ant + t_rcv )
              sefd        = (2*1380.00)/aot
           
           ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, za_deg_db , az_deg, za_deg )
           
           print("TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg] vs. limit = %.4f [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,ang_distance_deg,az_deg,za_deg,min_angular_distance_deg))
                        
           if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
               best_lst  = lst_db 
               best_aot  = aot 
               best_sefd = sefd
               best_id   = id
               min_angular_distance_deg = ang_distance_deg
           else :
               print("\t\tDB Record (azim,za) = (%.4f,%.4f) [deg] ignored due to angular distance too large ( %.4f deg )" % (azim_deg_db,za_deg_db,ang_distance_deg))
 
       if best_id > 0 :
          out_lst.append( best_lst )
          out_aot.append( best_aot )
          out_sefd.append( best_sefd ) 
          print("Use best matching record ID = %d (lst = %.4f [h], aot = %.4f , sefd = %.4f (minimum angular distance = %.4f [deg] , lst difference = %.4f [h])" % (best_id,best_lst,best_aot,best_sefd,min_angular_distance_deg,(best_lst-lst)))
       else :
          print("ERROR : no record found for lst = %.4f [h] and (ra,dec) = (%.4f,%.4f) [deg] -> (az,za) = (%.4f,%.4f) [deg]" % (lst_db,ra_deg, dec_deg, az_deg, za_deg))
 
    return (out_lst, out_aot, out_sefd)





# get sensitivity vs. time :
def get_sensitivity_timerange( az_deg , za_deg , freq_mhz, ux_start, ux_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=5.00, 
                             receiver_temperature=None,
                             ra_deg=None, dec_deg=None ) :
        
    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_uxtime_x, out_aot_x, out_sefd_x) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature, ra_deg=ra_deg, dec_deg=dec_deg )

    (out_uxtime_y, out_aot_y, out_sefd_y) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature, ra_deg=ra_deg, dec_deg=dec_deg )

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_uxtime_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_uxtime_x, out_sefd_y, out_uxtime_y )

    return ( numpy.array(out_uxtime_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_uxtime_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_uxtime_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )

def get_sensitivity_lstrange( az_deg , za_deg , freq_mhz, lst_start, lst_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None ) :

    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_lst_x, out_aot_x, out_sefd_x) = get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    (out_lst_y, out_aot_y, out_sefd_y) = get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_lst_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_lst_x, out_sefd_y, out_lst_y )


    return ( numpy.array(out_lst_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_lst_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_lst_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )

def get_sensitivity_radec_lstrange( ra_deg , dec_deg , freq_mhz, lst_start, lst_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None ) :

    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_lst_x, out_aot_x, out_sefd_x) = get_sensitivity_radec_lstrange_single_pol( ra_deg , dec_deg , freq_mhz, lst_start, lst_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    (out_lst_y, out_aot_y, out_sefd_y) = get_sensitivity_radec_lstrange_single_pol( ra_deg , dec_deg , freq_mhz, lst_start, lst_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_lst_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_lst_x, out_sefd_y, out_lst_y )


    return ( numpy.array(out_lst_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_lst_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_lst_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )



# get sensitivity map for a given LST and freq_mhz
def get_sensitivity_map( freq_mhz, lst_hours, 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_freq_resolution_mhz=10.00, output_file_base=None,
                             out_fitsname_base="sensitivity_map_", output_dir="./", receiver_temperature=None ) :
                       
    global debug_level                        
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )

    if conn is None :
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)

  
    # select closest LST :
    cur = conn.cursor()
    szSQL = "SELECT MIN(ABS(lst-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<=%.4f AND ABS(frequency_mhz-%.4f)<=%.4f" %  (lst_hours,lst_hours,db_lst_resolution, freq_mhz, db_freq_resolution_mhz )
    print("DEBUG SQL_best_lst : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_lst_distance = None
    for row in rows:
       if row[0] is not None :
          min_lst_distance = float( row[0] )
       else :
          print("ERROR : not data close to frequency = %.2f MHz and lst = %.2f h in the database" % (freq_mhz,lst_hours))
    
    if min_lst_distance is None :
       print("ERROR no records in the database exist closer than %.4f hours in LST at frequency %.2f MHz" % (db_lst_resolution,freq_mhz))
       return (None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)

    print("Best LST in database is closer than %.8f [hours]" % (min_lst_distance))
    min_lst_distance = min_lst_distance + 0.01

    # select closest FREQUENCY :
    cur = conn.cursor()
    szSQL = "SELECT MIN(ABS(frequency_mhz-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<=%.4f AND ABS(frequency_mhz-%.4f)<=%.4f" %  (freq_mhz, lst_hours, min_lst_distance, freq_mhz, db_freq_resolution_mhz )
    print("DEBUG SQL_best_freq : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_freq_distance = None
    for row in rows:
       print(row)
       min_freq_distance = float( row[0] )
    print("Best FREQ in database is closer than %.8f [MHz]" % (min_freq_distance))
    
    if min_freq_distance is None :
       print("ERROR no records in the database exist closer than %.4f MHz in FREQ at LST around %.4f [hours]" % (db_freq_resolution_mhz,lst_hours))
       return (None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
 
    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<=%.8f AND ABS(frequency_mhz-%.4f)<=%.8f ORDER BY za_deg, azim_deg ASC" %  (lst_hours,(min_lst_distance+0.01), freq_mhz, (min_freq_distance+0.1) )
    print("DEBUG SQL_main : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
 
    out_azim_x = []
    out_za_x   = []
    out_aot_x  = []
    out_sefd_x = []
    
    out_azim_y = []
    out_za_y   = []
    out_aot_y  = []
    out_sefd_y = []

    out_txt_x_f = None
    out_txt_y_f = None
    out_txt_i_f = None

    out_txt_filename_X = None
    out_txt_filename_Y = None
    out_txt_filename_I = None
    if output_file_base is not None :
       out_txt_filename_X = output_dir + "/" + out_fitsname_base + "_X.txt"
       out_txt_filename_Y = output_dir + "/" + out_fitsname_base + "_Y.txt"
       out_txt_filename_I = output_dir + "/" + out_fitsname_base + "_I.txt"
       
       out_txt_x_f = open( out_txt_filename_X , "w" )
       out_txt_y_f = open( out_txt_filename_Y , "w" )
       out_txt_i_f = open( out_txt_filename_I , "w" )
       
       print("DEBUG : saving output text files to %s and %s" % (out_txt_filename_X,out_txt_filename_Y))
       
       out_txt_x_f.write( code_version_line )
       out_txt_y_f.write( code_version_line )
       out_txt_i_f.write( code_version_line )
       
       line = "# AZIM[deg] ZA[deg] A/T[m^2/K]"
       out_txt_x_f.write( line + "    for X polarisation\n" )
       out_txt_y_f.write( line + "    for Y polarisation\n" )
       out_txt_i_f.write( line + "    for Stokes I polarisation\n" )              
        
    for row in rows:
        if debug_level >= 2 :
           print(row)
        
        id = int( row[0] )
        azim_deg_db = float( row[1] )
        za_deg_db   = float( row[2] )
        freq_mhz    = float( row[3] )
        pol         = row[4]
        lst_db      = float( row[5] )
        unixtime    = float( row[6] )
        gpstime     = float( row[7] )
        aot         = float( row[8] )
        sefd        = numpy.NaN
        if aot > 0 :
            sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        if receiver_temperature is not None and receiver_temperature >= 0.00 :
           t_rcv = receiver_temperature
           aot = a_eff / ( t_ant + t_rcv )
           sefd        = (2*1380.00)/aot
        
        if debug_level >= 2 :
           print("TEST : %d , freq_mhz = %.4f [MHz]" % (id,freq_mhz))
        
        line = "%.4f %.4f %.8f\n" % (azim_deg_db,za_deg_db,aot)
        
        if pol == "X" :
           out_azim_x.append( azim_deg_db )
           out_za_x.append( za_deg_db )
           out_aot_x.append( aot )
           out_sefd_x.append( sefd )
           
           if out_txt_x_f is not None : 
               out_txt_x_f.write( line )                      
        elif pol == "Y" :
           out_azim_y.append( azim_deg_db )
           out_za_y.append( za_deg_db )
           out_aot_y.append( aot )
           out_sefd_y.append( sefd )            
           
           if out_txt_y_f is not None : 
               out_txt_y_f.write( line )                      
        else :
           print("ERROR : unknown polarisation = %s" % (pol))
 
    ( out_freq_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_azim_x, out_sefd_y, out_za_x, out_file=out_txt_i_f  )   
    
    
    if out_txt_x_f is not None :
       out_txt_x_f.close()

    if out_txt_y_f is not None :
       out_txt_y_f.close()

    if out_txt_i_f is not None :
       out_txt_i_f.close()


    return ( numpy.array(out_azim_x), numpy.array(out_za_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_azim_y), numpy.array(out_za_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_azim_y), numpy.array(out_za_y), numpy.array(out_aot_i) , numpy.array(out_sefd_i),
             out_txt_filename_X, out_txt_filename_Y, out_txt_filename_I )



# 
def get_sensitivity_azzalstrange( az_deg , za_deg , freq_mhz, lst_start_h, lst_end_h , 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, freq_resolution_mhz=5.00,
                             receiver_temperature=None ) :
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
   
    if conn is None :    
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None)



    # get requested data :
    cur = conn.cursor()
    # Condition %.4f<0.1 means that if za_deg is close to 0 (zenith) azimuth does not matter !
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE lst between %.4f and %.4f AND ABS(za_deg-%.4f)<=%.4f AND (ABS(azim_deg-%.4f)<=%.4f OR %.4f<0.1) AND ABS(frequency_mhz-%.4f)<=%.4f ORDER BY lst ASC" %  (lst_start_h,lst_end_h,za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, za_deg, freq_mhz, freq_resolution_mhz )
    cur.execute( szSQL )
    rows = cur.fetchall()
 
    out_lst_x = []
    out_aot_x  = []
    out_sefd_x = []
    
    out_lst_y = []
    out_aot_y  = []
    out_sefd_y = []
        
    for row in rows:
        print(row)
        
        id = int( row[0] )
        azim_deg_db = float( row[1] )
        za_deg_db   = float( row[2] )
        freq_mhz    = float( row[3] )
        pol         = row[4]
        lst_db      = float( row[5] )
        unixtime    = float( row[6] )
        gpstime     = float( row[7] )
        aot         = float( row[8] )
        sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        if receiver_temperature is not None and receiver_temperature >= 0.00 :
           t_rcv = receiver_temperature
           aot = a_eff / ( t_ant + t_rcv )
           sefd        = (2*1380.00)/aot
        
        print("TEST : %d , freq_mhz = %.4f [MHz]" % (id,freq_mhz))
        
        if pol == "X" :
           out_lst_x.append( lst_db )
           out_aot_x.append( aot )
           out_sefd_x.append( sefd ) 
        elif pol == "Y" :
           out_lst_y.append( lst_db )
           out_aot_y.append( aot )
           out_sefd_y.append( sefd )            
        else :
           print("ERROR : unknown polarisation = %s" % (pol))
 

    return ( numpy.array(out_lst_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_lst_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )



def plot_sensitivity_vs_time( uxtime_x, aot_x, uxtime_y, aot_y,  unixtime_start, unixtime_end, azim_deg, za_deg, freq_mhz,
                              output_file_base=None, point_x='go', point_y='rx',
                              min_ylimit=0.00, max_ylimit=2.00,
                              uxtime_i=None, aot_i=None , point_i='b+' ,
                              fig_size_x=20, fig_size_y=10, info=None, save_output_path="./", save_text_file=False, out_dir="./" ) :
                     
   global plot_requirements
   global plot_braun2019

   legend_list = []

#   legend_location = "upper center" # "upper right"
   legend_location = "upper right"
   
   # MAX A/T :
   max_aot = 0.00    
   if aot_x is not None :
      if max(aot_x) > max_aot :
         max_aot = max(aot_x)

   if aot_y is not None :
      if max(aot_y) > max_aot :
         max_aot = max(aot_y)

   if aot_i is not None :
      if max(aot_i) > max_aot :
         max_aot = max(aot_i)


   # conversion of unix time to UTC :                               
   x_utc = []
   y_utc = []
   for ux in uxtime_x :
#      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      utc = datetime.utcfromtimestamp(ux)
      print("DEBUG : %d -> %s" % (ux,utc))
      x_utc.append( utc )
            
   for ux in uxtime_y :
#      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      utc = datetime.utcfromtimestamp(ux)
      y_utc.append( utc )
   
   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   fig = plt.gcf()
#   fig.axes[0].xaxis.set_major_formatter( md.DateFormatter('%Y') )
   
   if uxtime_x is not None and aot_x is not None :
      ax_x =  plt.plot( x_utc, aot_x, point_x )
      legend_list.append( 'X polarisation' )

      # ax_x.tick_params(which='major', length=10, width=20, direction='inout')
#      ax_x.set_ylim(( min_ylimit,  max_ylimit ))
   
   if uxtime_y is not None and aot_y is not None :
      ax_y = plt.plot( y_utc, aot_y, point_y )
      legend_list.append( 'Y polarisation' )

   if uxtime_i is not None and aot_i is not None :
      ax_i = plt.plot( x_utc, aot_i, point_i )
      legend_list.append( 'Stokes I' )
      
      max_i = max(aot_i)
      if max_i > max_ylimit :
         max_ylimit = max_i*1.1 # max_i + 10% max_i
#      ax_y.set_ylim(( min_ylimit,  max_ylimit ))
 
   if plot_requirements :
      sens_req_zenith = []
      sens_req_avg45 = []

      for ux in uxtime_x :
         print("DEBUG : getting requirement for uxtime = %.2f" % (ux))
         sens_zenith = lfaa_requirements.lfaa_per_station( freq_mhz, subversion="zenith" )
         sens_avg45  =  lfaa_requirements.lfaa_per_station( freq_mhz, subversion="AvgAboveElev45deg" )
         print("DEBUG : getting requirement DONE")
         
         sens_req_zenith.append( sens_zenith )
         sens_req_avg45.append( sens_avg45 )
         
         print("DEBUG : appended")
         
      sens_req_zenith = numpy.array( sens_req_zenith )
      sens_req_avg45 = numpy.array( sens_req_avg45 )

      # 
      print("DEBUG : plotting SKA-Low requrirements at ZENITH (%d data points)..." % (len(uxtime_x)))
      # sens_req_zenith,
      plt.plot( x_utc, sens_req_zenith, linestyle='dashed', color='orange', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : plotting SKA-Low requrirements at ZENITH DONE")
      legend_list.append('SKA-Low requirement at zenith')
      
      plt.plot( x_utc, sens_req_avg45, linestyle='dotted', color='red', linewidth=2, markersize=12 ) # marker='o'
      legend_list.append('SKA-Low requirement (average above $45^o$)')
 
      # plt.legend(('X polarisation','Y polarisation','Stokes I'), loc='upper right' , fontsize=20) 

   if plot_braun2019 :
      sens_braun = []
 
      for utc in x_utc :
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019) for utc = %s" % (utc))
         sens_value = sensitivity_braun2019.lfaa_per_station( freq_mhz )
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019)  DONE , A/T = %.4f [m^2/K]" % (sens_value))

         sens_braun.append( sens_value )
         print("DEBUG : appended")

      sens_braun = numpy.array( sens_braun )

      # 
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun ...")
      plt.plot( x_utc, sens_braun, linestyle='dashed', color='magenta', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun DONE")
      legend_list.append('Table 9 in Braun et al. (2019)')


   # legend :
   if uxtime_x is not None and aot_x is not None and uxtime_y is not None and aot_y is not None :
      if uxtime_i is not None and aot_i is not None :
         plt.legend( legend_list, loc=legend_location, fontsize=20)
      else :
         plt.legend( leged_list, loc=legend_location,  fontsize=20 )
   else :
      if uxtime_x is not None and aot_x is not None :
         plt.legend( legend_list, loc=legend_location,  fontsize=20 )
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=legend_location,handles=[uxtime_y, uxtime_y])
   
   plt.xlabel('Time [UTC]' , fontsize=30 )
   plt.ylabel('Sensitivity A / T [m$^2$/K]' , fontsize=30 )

   ax_list = fig.axes
   plt.gcf().autofmt_xdate()
   
   ticks = plt.xticks
   x_start, x_end = ax_list[0].get_xlim()
   print("DEBUG : x_start = %.8f , x_end = %.8f" % (x_start,x_end))
   print("DEBUG : x_end - x_start= %.2f" % (x_end - x_start))
   
   # change format of the date in X-axis :
   xaxis_fmt = md.DateFormatter('%d/%m %H:%M')
   ax_list[0].xaxis.set_major_formatter(xaxis_fmt)
#   ax_list[0].xaxis.set_major_locator(md.MinuteLocator(interval=60))
   ax_list[0].xaxis.set_major_locator(md.HourLocator(interval=4))   

#   if ( x_end - x_start ) >= 20.00 :
#      ax_list[0].xaxis.set_ticks(numpy.arange( x_start, x_end, 8)) # 4 -> 8
#   else :
#      ax_list[0].xaxis.set_ticks(numpy.arange( x_start, x_end, 2))
   
   
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M')
#   ax_list[0].xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))
   max_ylimit0 = max_ylimit
   max_ylimit = max_ylimit*1.1 # + 10%
   
   if max_aot > (max_ylimit-1.8) :
      max_ylimit = max_aot + 1.8
   
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
   if info is not None :
      # place a text box in upper left in axes coords
      min_x = min( x_utc )
      text( min_x, max_ylimit*0.85, info , fontsize=30 ) # , transform=ax.transAxes, verticalalignment='top', bbox=props)
   
   plt.grid()
   
   if output_file_base is not None :
      outfile = ( "%s/%s.png" % (save_output_path,output_file_base) )
      plt.savefig( outfile )
   
   plt.show()
   
      # save txt file :
   out_txt_filename = None
   if save_text_file : # saving of full (interpolated) map as in FITS file is turned of as the files are too large :
      out_txt_filename = out_dir + "/" + output_file_base + "_X.txt"
      out_txt_f = open( out_txt_filename , "w" )
      out_txt_f.write( code_version_line )
      
      line = "# UXTIME A/T_X[m^2/K] A/T_Y[m^2/K] A/T_I[m^2/K]\n"
      out_txt_f.write( line )
      
      if len(uxtime_x) == len(uxtime_y) :      
         for i in range(0,len(uxtime_x)) :
            (aot_i,sefd_i) = aotxy2aoti( aot_x[i], aot_y[i] )
            line = "%.2f %.8f %.8f %.8f\n" % (uxtime_x[i],aot_x[i],aot_y[i],aot_i)
            out_txt_f.write( line )
         
      out_txt_f.close()
 
      print("Sensitivity map saved to text file %s" % (out_txt_filename))

   
def plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y,  lst_start, lst_end, azim_deg, za_deg, freq_mhz,
                              output_file_base=None, point_x='go', point_y='rx',
                              min_ylimit=0.00, max_ylimit=2.00,
                              do_show = True, save_output_path="./",
                              lst_i=None, aot_i=None, point_i='b+',
                              fig_size_x=20, fig_size_y=10, info=None,
                              ra_deg=None, dec_deg=None ) :

   global web_interface_initialised
   global plot_requirements
   global plot_braun2019

   legend_list = []

   max_aot = 0.00   
   min_lst = 0.00
   
   if info is None :
      if ra_deg is not None and dec_deg is not None :
         info = "Frequency %.2f MHz, (RA,Dec) = (%.2f$^o$, %.2f$^o$)" % (freq_mhz,ra_deg,dec_deg)
      else :
         info = "Frequency %.2f MHz, Azimuth = %.2f$^o$, ZA = %.2f$^o$" % (freq_mhz,azim_deg,za_deg)
   
   if lst_x is not None :
      min_lst = min(lst_x)
   
   if aot_x is not None :
      if max(aot_x) > max_aot :
         max_aot = max(aot_x)
         
   if aot_y is not None :
      if max(aot_y) > max_aot :
         max_aot = max(aot_y)

   if aot_i is not None :
      if max(aot_i) > max_aot :
         max_aot = max(aot_i)
   
   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   if lst_x is not None and aot_x is not None :
      ax_x =  plt.plot( lst_x, aot_x, point_x )
      legend_list.append( 'X polarisation' )
   
   if lst_y is not None and aot_y is not None :
      ax_y = plt.plot( lst_y, aot_y, point_y )
      legend_list.append( 'Y polarisation' )

   if lst_i is not None and aot_i is not None :
      ax_i = plt.plot( lst_i, aot_i, point_i )
      
      max_i = max( aot_i )
      legend_list.append( 'Stokes I' )
#      if max_i > max_ylimit :
#         max_ylimit = max_i * 1.1 # max_i + 10% of max_i


   if plot_requirements :
      lst_req = []
      sens_req_zenith = []
      sens_req_avg45 = []
      
      for lst in lst_x :
         print("DEBUG : getting requirement for lst = %.4f [hours]" % (lst))
         sens_zenith = lfaa_requirements.lfaa_per_station( freq_mhz, subversion="zenith" )
         sens_avg45  =  lfaa_requirements.lfaa_per_station( freq_mhz, subversion="AvgAboveElev45deg" )
         print("DEBUG : getting requirement DONE")
         
         lst_req.append( lst )
         sens_req_zenith.append( sens_zenith )
         sens_req_avg45.append( sens_avg45 )
         
         print("DEBUG : appended")
         
      lst_req = numpy.array( lst_req )   
      sens_req_zenith = numpy.array( sens_req_zenith )
      sens_req_avg45 = numpy.array( sens_req_avg45 )

      # 
      print("DEBUG : plotting SKA-Low requrirements at ZENITH ...")         
      plt.plot( lst_req, sens_req_zenith, linestyle='dashed', color='orange', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : plotting SKA-Low requrirements at ZENITH DONE")
      legend_list.append('SKA-Low requirement at zenith')
      
      plt.plot( lst_req, sens_req_avg45, linestyle='dotted', color='red', linewidth=2, markersize=12 ) # marker='o'
      legend_list.append('SKA-Low requirement (average above $45^o$)')
      
      # plt.legend(('X polarisation','Y polarisation','Stokes I'), loc='upper right' , fontsize=20) 

   if plot_braun2019 :
      lst_braun = []
      sens_braun = []
      
      for lst in lst_x :
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019) for lst = %.4f [hours]" % (lst))
         sens_value = sensitivity_braun2019.lfaa_per_station( freq_mhz )
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019)  DONE , A/T = %.4f [m^2/K]" % (sens_value))

         lst_braun.append( lst )
         sens_braun.append( sens_value )
         print("DEBUG : appended")

      lst_braun = numpy.array( lst_braun )   
      sens_braun = numpy.array( sens_braun )

      # 
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun ...")
      plt.plot( lst_braun, sens_braun, linestyle='dashed', color='magenta', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun DONE")
      legend_list.append('Table 9 in Braun et al. (2019)')


   # legend :
   if lst_x is not None and aot_x is not None and lst_y is not None and aot_y is not None :
      if lst_i is None and aot_i is None :      
         plt.legend( legend_list , loc='upper right' , fontsize=20 )
      else :
         plt.legend( legend_list, loc='upper right' , fontsize=20 )
   else :
      if lst_x is not None and aot_x is not None :
         plt.legend( legend_list , loc='upper right' , fontsize=20 )
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[lst_y, lst_y])
   
   plt.xlabel('Local sidereal time [hours]' , fontsize=30 )
   plt.ylabel('Sensitivity A / T [m$^2$/K]' , fontsize=30 )
   plt.gcf().autofmt_xdate()
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M:%S')
#   ax.xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))

   if max_aot > (max_ylimit-0.5) :
      max_ylimit = max_aot + 0.5
      
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
   if info is not None :
      # place a text box in upper left in axes coords
      text( min_lst*0.75, max_ylimit*1.05, info , fontsize=25 ) # , transform=ax.transAxes, verticalalignment='top', bbox=props)
   
   plt.grid()
   
   png_image_path = None
   if output_file_base is not None :
      png_image_path = ( "%s/%s.png" % ( save_output_path , output_file_base) )
      
      plt.savefig( png_image_path )
      print("Saved output image to file %s" % (png_image_path))
      
      if web_interface_initialised :
#         sio = StringIO() # use io.StringIO() in Python 3x
#         sio = io.StringIO()
         buf = BytesIO()
         plt.savefig( buf, format="png")
#         buf.seek(0)
#         string = base64.b64encode(buf.read())
         
#         print("web_interface_initialised = True -> returning image = %s" % (buf))

         # WORKS : see https://stackoverflow.com/questions/52368870/display-matplotlib-image-on-html-page-using-django         
         return (png_image_path,buf)

   if do_show :   
      plt.show()
   

   return (png_image_path)

def plot_sensitivity( freq_x, aot_x, freq_y, aot_y, output_file_base=None, point_x='go', point_y='rx', freq_i=None, aot_i=None, point_i='b+', min_ylimit=0.00, max_ylimit=2.00, info=None,
                      fig_size_x=20, fig_size_y=10, save_output_path="./" ):

   global web_interface_initialised
   global plot_requirements
   global plot_braun2019
   
   print("DEBUG : inside plot_sensitivity , plot_requirements = %s" % (plot_requirements))
   
   max_aot = 0.00
   min_freq = min(freq_x)
   
   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   ax = None
   legend_list = []
   if freq_x is not None and aot_x is not None :
      plt.plot( freq_x, aot_x, point_x )
      legend_list.append( 'X polarisation' )
      
      if max(aot_x) > max_aot :
         max_aot = max(aot_x)
   
   if freq_y is not None and aot_y is not None :
      plt.plot( freq_y, aot_y, point_y )
      legend_list.append( 'Y polarisation' )
      
      if max(aot_y) > max_aot :
         max_aot = max(aot_y)


   if freq_i is not None and aot_i is not None :
      plt.plot( freq_i, aot_i, point_i )
      legend_list.append( 'Stokes I' )

      if max(aot_i) > max_aot :
         max_aot = max(aot_i)
         
   if plot_requirements :
      freq_req = []
      sens_req_zenith = []
      sens_req_avg45 = []
      
      for freq_mhz in freq_x :
         print("DEBUG : getting requirement for freq = %.4f MHz" % (freq_mhz))
         sens_zenith = lfaa_requirements.lfaa_per_station( freq_mhz, subversion="zenith" )
         sens_avg45  =  lfaa_requirements.lfaa_per_station( freq_mhz, subversion="AvgAboveElev45deg" )
         print("DEBUG : getting requirement DONE")
         
         freq_req.append( freq_mhz )
         sens_req_zenith.append( sens_zenith )
         sens_req_avg45.append( sens_avg45 )
         
         print("DEBUG : appended")
         
      freq_req = numpy.array( freq_req )   
      sens_req_zenith = numpy.array( sens_req_zenith )
      sens_req_avg45 = numpy.array( sens_req_avg45 )

      # 
      print("DEBUG : plotting SKA-Low requrirements at ZENITH ...")         
      plt.plot( freq_req, sens_req_zenith, linestyle='dashed', color='orange', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : plotting SKA-Low requrirements at ZENITH DONE")
      legend_list.append('SKA-Low requirement at zenith')
      
      plt.plot( freq_req, sens_req_avg45, linestyle='dotted', color='red', linewidth=2, markersize=12 ) # marker='o'
      legend_list.append('SKA-Low requirement (average above $45^o$)')
      
      # plt.legend(('X polarisation','Y polarisation','Stokes I'), loc='upper right' , fontsize=20) 

   if plot_braun2019 :
      freq_braun = []
      sens_braun = []
      
      for freq_mhz in freq_x :
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019) for freq = %.4f MHz" % (freq_mhz))
         sens_value = sensitivity_braun2019.lfaa_per_station( freq_mhz )
         print("DEBUG : getting Table 9 (avg above 45deg) from Braun et al. (2019)  DONE , A/T = %.4f [m^2/K]" % (sens_value))

         freq_braun.append( freq_mhz )
         sens_braun.append( sens_value )
         print("DEBUG : appended")

      freq_braun = numpy.array( freq_braun )   
      sens_braun = numpy.array( sens_braun )

      # 
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun ...")
      plt.plot( freq_braun, sens_braun, linestyle='dashed', color='magenta', linewidth=2, markersize=12 ) # marker='o'
      print("DEBUG : over-plotting Table 9 (avg above 45deg) from Braun DONE")
      legend_list.append('Table 9 in Braun et al. (2019)')


   # legend :
   if freq_x is not None and aot_x is not None and freq_y is not None and aot_y is not None :
      if freq_i is not None and aot_i is not None :
         plt.legend(legend_list, loc='upper right' , fontsize=20)    
      else :
         plt.legend(legend_list, loc='upper right' , fontsize=20)
   else :
      if freq_x is not None and aot_x is not None :
         plt.legend(legend_list, loc='upper right' , fontsize=20)
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[freq_y, aot_y],fontsize=20)
   
   plt.xlabel('Frequency (MHz)' , fontsize=30 )
   plt.ylabel('Sensitivity A/T [m$^2$/K]' , fontsize=30 )
   plt.grid()

#   if max_aot < 0.5 :
#      max_ylimit = 0.5
      
   max_ylimit0 = max_ylimit
   max_ylimit = max_ylimit*1.1 # + 10%

   # have at least 0.5 buffer at the top :
   if max_aot > ( max_ylimit - 0.8 ) :
      max_ylimit = max_aot + 0.8
   
   plt.ylim(( min_ylimit,  max_ylimit )) 
   
   if info is not None :
      # place a text box in upper left in axes coords
      text( min_freq*0.75, max_ylimit*0.95, info , fontsize=25 ) # , transform=ax.transAxes, verticalalignment='top', bbox=props)

   
   outfile = None
   if output_file_base is not None :
      outfile = ( "%s/%s.png" % (save_output_path,output_file_base) )
      print("Saving file %s ..." % (outfile))
      plt.savefig( outfile )

      print("Saved output image to file %s" % (outfile))
      
      if web_interface_initialised :
         buf = BytesIO()
         plt.savefig( buf, format="png")
         # WORKS : see https://stackoverflow.com/questions/52368870/display-matplotlib-image-on-html-page-using-django         
         return (outfile,buf)
   
   plt.show()
   
   return( outfile , None )

def get_radius_list( za_rad, za_list ) :
   size_x = za_rad.shape[0]
   size_y = za_rad.shape[1]
   
   radius_px_list = []
   for za_radius in za_list :
      min_diff = 1e10
      pixel_radius = -1
      for x in range(int(size_x/2),int(size_x)) :
         za_deg = za_rad[x,int(size_y/2)]*(180.00/math.pi)
         
         diff = math.fabs(za_deg - za_radius)
         if diff < min_diff :
            min_diff = diff
            pixel_radius = x - int(size_x/2)
            
      if pixel_radius > 0 :
         print("DEBUG : for za=%.4f [deg] radius in pixels = %d" % (za_radius,pixel_radius))
      else :
         print("ERROR : for za=%.4f [deg] radius in pixels is %d" % (za_radius,pixel_radius))
         
      radius_px_list.append( pixel_radius )
         
   return radius_px_list
      

def plot_sensitivity_map( azim_deg, za_deg, aot, out_fitsname_base="sensitivity" , save_text_file=False, do_plot=False, freq_mhz=0.00, lst_h=0.00, pol="Unknown", out_dir="./", s=3.5, vmax_value=2.00 ) :
   from scipy.interpolate import SmoothSphereBivariateSpline    
   global web_interface_initialised


   azim_rad = azim_deg*(numpy.pi/180.0)
   za_rad   = za_deg*(numpy.pi/180.0)
   lut = SmoothSphereBivariateSpline( za_rad, azim_rad, aot , s=s)
#   fine_lats = numpy.linspace(0., numpy.pi/2.00,256)
#   fine_lons = numpy.linspace(0., 2 * numpy.pi, 256)
#   data_smth = lut(fine_lats, fine_lons)
   
#   plt.figure()
#   plt.plot( freq_x, aot_x )
#   ax=fig.add_subplot(111)
#   ax.imshow(data_smth, interpolation='nearest')

#   def makeAZZA(npix=256,projection='SIN',azim_from_north=False,return_all=False):
   (az_rad,za_rad,x,y,z,d) = beam_tools.makeAZZA(npix=256,projection='SIN',azim_from_north=True, return_all=True )    
   
   
   sensitivity = numpy.zeros( az_rad.shape )
   max_sens = 0
   for i in range(0,za_rad.shape[0]) :
      for j in range(0,za_rad.shape[1]) :
         a_rad = az_rad[ i , j ]
         z_rad = za_rad[ i , j ]
         
         aot_value = lut( z_rad, a_rad )
         if aot_value >= 0 :
            # ignore non-physical negative values :
            sensitivity[ i , j ] = aot_value
            
            if aot_value > max_sens :
               max_sens = aot_value
         

   pol_string = "_%s" % (pol)
   sens_fits_file = out_dir + "/" + out_fitsname_base + pol_string + ".fits"
   azim_fits_file = out_dir + "/" + out_fitsname_base + pol_string + "_azim_deg.fits"
   za_fits_file   = out_dir + "/" + out_fitsname_base + pol_string + "_za_deg.fits"
   
   fits_beam.save_fits( sensitivity , sens_fits_file )
   fits_beam.save_fits( az_rad*(180.00/numpy.pi) , azim_fits_file )
   fits_beam.save_fits( za_rad*(180.00/numpy.pi) , za_fits_file )
         
   # save txt file :
   out_txt_filename = None
   if save_text_file : # saving of full (interpolated) map as in FITS file is turned of as the files are too large :
      out_txt_filename = out_dir + "/" + out_fitsname_base + pol_string + ".txt"
      out_txt_f = open( out_txt_filename , "w" )
      out_txt_f.write( code_version_line )
      line = "# AZIM[deg] ZA[deg] A/T[m^2/K]\n"
      out_txt_f.write( line )
      for i in range(0,za_rad.shape[0]) :
         for j in range(0,za_rad.shape[1]) :
            a_rad = az_rad[ i , j ]
            z_rad = za_rad[ i , j ]
 
            aot_value = lut( z_rad, a_rad )
            if aot_value >= 0 :
               line = "%.4f %.4f %.8f\n" % (a_rad*(180.00/math.pi),z_rad*(180.00/math.pi),aot_value)
               out_txt_f.write( line )
     
      out_txt_f.close()
      
      print("Sensitivity map saved to text file %s" % (out_txt_filename))

   pngfile = None
   buf     = None
   
   if do_plot :
      # see /home/msok/ska/aavs/aavs0.5/trunk/simulations/FEKO/beam_models/MWA_EE/MWAtools_pb/primarybeammap_local.py
      figsize=16
      fig=plt.figure(figsize=(figsize,0.6*figsize),dpi=300)
      axis('on')
      ax1=fig.add_subplot(1,1,1,polar=False)
    
      axis('off')
      ax2=fig.add_subplot(1,1,1,polar=True, frameon=False)
      ax2.set_theta_zero_location("N")
      ax2.set_theta_direction(-1)
      ax2.patch.set_alpha(0.0) # 0.0 -> 1.0?
      ax2.tick_params(color='0.5', labelcolor='0.5', pad=-40 ) # color='0.5', labelcolor='0.5'
      for spine in ax2.spines.values():
          spine.set_edgecolor('0.5')
          
      # WARNING : 2021-11-15 :
      # Grid is not really SIN projection -> turned off for now            
#      ax2.grid(which='major', color='0.5') 
      ax2.grid(b=None)

      # hide distance from the center
      # WARNING : 2021-11-15 : Grid is not really SIN projection -> turned off for now 
      for xlabel_i in ax2.get_yticklabels():
         xlabel_i.set_visible(False)
    
      #Beamsky example:
#      vmax = 1.75
#      vmax = max_sens
      vmax=vmax_value
      im=ax1.imshow( sensitivity , interpolation='none', vmax=vmax )

      # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Circle.html      
      # add circles with za=15, 30, 45, 75 degrees :
      za_circle_list = (15,30,45,75)
      radius_list = get_radius_list( za_rad, za_circle_list )
#      radius_list = get_radius_list( za_rad, (10,20,30,40,50) )
      circle_index = 0
      circle_color="red" # was lightgrey
      for radius in radius_list :
         circ = Circle((az_rad.shape[0]/2,az_rad.shape[1]/2),radius,color=circle_color,fill=False,linewidth=0.3,linestyle="--")
         ax1.add_patch(circ)
         
         str = "$%d^o$" % (za_circle_list[circle_index])
#         plt.text( az_rad.shape[0]/2,az_rad.shape[1]/2 + za_circle_list[circle_index], str )
#         ax1.text( az_rad.shape[0]/2 + radius_list[circle_index] - 10, az_rad.shape[1]/2 - (circle_index+1)*5 , str, color=circle_color )

         angle_deg=35.00         
         x_offset = radius_list[circle_index]*math.cos( angle_deg*math.pi/180.00 )
         y_offset = radius_list[circle_index]*math.sin( angle_deg*math.pi/180.00 )

         ax1.text( az_rad.shape[0]/2 + x_offset, az_rad.shape[1]/2 - y_offset , str, color=circle_color, fontsize=30 )
         circle_index += 1

      # TEST :          
#      circ = Circle((az_rad.shape[0]/2,az_rad.shape[1]/2),120,color="red",fill=False,linewidth=0.3)
#      ax1.add_patch(circ)   

      #Add colorbar on own axis
      cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
      cbar=fig.colorbar(im, cax=cbar_ax)
      cbar.set_label("A / T [m$^2$/K]", fontsize=30, labelpad=10)
#      cbar.ax.set_ylim( 0 , 1.75 )
      
      # ticks :
      # ax1.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
      # ax1.set_rticks([])
 
      title = "Sensitivity map at frequency = %.2f MHz at lst = %.2f [h] (%s)" % (freq_mhz,lst_h,pol)
      ax1.set_title( title , fontsize=17, pad=15 )
      pol_str = "_%s" % (pol)
      pngfile = out_dir + "/" + out_fitsname_base + pol_str + ".png"
      fig.savefig( pngfile )
      print("Saved file %s" % (pngfile))
      plt.show()
      
      
      if web_interface_initialised :
         buf = BytesIO()
         plt.savefig( buf, format="png")
         # WORKS : see https://stackoverflow.com/questions/52368870/display-matplotlib-image-on-html-page-using-django         

   
   return (lut,sensitivity,sens_fits_file,azim_fits_file,za_fits_file,out_txt_filename,pngfile,buf)
 

 
def parse_options(idx):
   global debug_level

   usage="Usage: %prog [options]\n"
   usage += "Different plots / data can be obtainted:\n"
   usage += " 1/ to plot AoT vs. freq. for a given pointing direction (az,za) and lst time use options --azim_deg, --za_deg, --lst_hours\n"
   usage += " 2/ to plot AoT vs. time for a specified frequency use options --freq_mhz, --lst_start, --lst_end\n"
   usage += " 3/ to plot sensitivity map over the entire sky at a specified LST and frequency specify options : --lst_hours and --freq_mhz\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # General parameters :
   parser.add_option('-s','--station_name','--station',dest="station_name",default="AAVS2", help="Station name [default %default]")

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-p','--plot','--do_plot','--do_plots',action="store_true",dest="do_plot",default=False, help="Plot")
   parser.add_option('-a','--azim','--azim_deg',dest="azim_deg",default=None, help="Pointing direction azimuth in degrees [default %default]",metavar="float",type="float")
   parser.add_option('-z','--za','--za_deg',dest="za_deg",default=None, help="Pointing direction zenith distance in degrees [default %default]",metavar="float",type="float")
   parser.add_option('--pointing_ra_deg','--ra','--ra_deg',dest='pointing_ra_deg',default=None,help='Pointing RA [deg]',type=float)
   parser.add_option('--pointing_dec_deg','--dec','--dec_deg',dest='pointing_dec_deg',default=None,help='Pointing DEC [deg]',type=float)
   parser.add_option('--object','--object_name',dest='object_name',default=None,help='Object name [default not specified')
   parser.add_option('-l','--lst','--lst_hours',dest="lst_hours",default=None, help="Local sidereal time in hours [default %default]",metavar="float",type="float")

   # specific frequency in MHz 
   parser.add_option('-f','--freq_mhz','--frequency_mhz',dest="freq_mhz",default=None, help="Specific frequency in MHz [default %default]",metavar="float",type="float")

   # specify lst range :
   parser.add_option('--lst_start','--lst_start_hours',dest="lst_start_hours",default=None, help="Start time in Local sidereal time in hours [default %default]",metavar="float",type="float")
   parser.add_option('--lst_end','--lst_end_hours',dest="lst_end_hours",default=None, help="End time in Local sidereal time in hours [default %default]",metavar="float",type="float")

   # specify Unix time range :
   parser.add_option('--ux_start','--unixtime_start',dest="unixtime_start",default=None, help="Start time in unixtime [default %default]",metavar="float",type="float")
   parser.add_option('--ux_end','--unixtime_end',dest="unixtime_end",default=None, help="End time in unixtime [default %default]",metavar="float",type="float")
   parser.add_option('--ux_interval','--unixtime_interval','--interval',dest="unixtime_interval",default=86400, help="Interval in unixtime [default %default]",metavar="float",type="float")   

   # specify UT range :
   parser.add_option('--ut_start','--utc_start','--start_utc',dest="ut_start",default=None, help="Start time in UTC [default %default]")
   parser.add_option('--ut_end','--utc_end','--end_utc',dest="ut_end",default=None, help="End time in UTC [default %default]")

   # time step for ranges :   
   parser.add_option('--timestep','--step_seconds','--step_sec',dest="step_seconds",default=300, help="Time step in seconds for both unix/utc-string time ranges [default %default]",metavar="float",type="float")

   # output file :
   parser.add_option('-o','--out_file','--outfile','--outout_file',dest="output_file",default=None, help="Full path to output text file basename (X or Y is added at the end) [default %default]" )
   parser.add_option('--save_text','--save_text_file','--text',action="store_true",dest="save_text_file",default=False, help="Save text file - for some options is always enabled, but for all-sky sensitivity map it's not [default %]")   
   
 
   # providing external values of T_receiver :
   parser.add_option('--receiver_temperature','--receiver_temp','--t_rcv','--T_rcv',dest="receiver_temperature",default=None, help="Externally provided receiver temperature (overwrites the value in the database) [default %default]",metavar="float",type="float")
   parser.add_option('--trcv_file','--receiver_temp_file',dest="receiver_temp_file",default=None, help="Text file with receiver temperature vs. frequency [default %default]")
   
   # debug / verbose level :
   parser.add_option('--debug_level','--verb_level','--verb','--debug',dest="debug_level",default=0, help="Debug / verbosity level [default %default]",type="int")
      

   (options, args) = parser.parse_args(sys.argv[idx:])
   
   debug_level = options.debug_level
   
   
   # unix time range has priority over UTC string range :
   # If unix time range is not filled in, but UTC string range is -> use unixtime range to later only refer in the code to UNIXTIME ranges 
   if options.unixtime_start is None or options.unixtime_end is None :
      if options.ut_start is not None and options.ut_end is not None :
         t_start_utc = Time( options.ut_start, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))
         t_end_utc   = Time( options.ut_end, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))
      
         if options.unixtime_start is None :
            options.unixtime_start = t_start_utc.unix

         if options.unixtime_end is None :
            options.unixtime_end = t_end_utc.unix

   if options.unixtime_start is not None and (options.unixtime_end is None and options.unixtime_interval is not None) :
      options.unixtime_end = options.unixtime_start + options.unixtime_interval
      
   if options.unixtime_start is not None and  options.unixtime_end is not None :
      if options.unixtime_interval is None :
         options.unixtime_interval = ( options.unixtime_end - options.unixtime_start )

   if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
      uxtime = options.unixtime_start
      if uxtime is None :
         if options.lst_hours is None and (options.lst_start_hours is None or options.lst_end_hours is None) :
            print("ERROR : not unixtime nor LST not specified use option --unixtime_start or --uxtime_start or --lst or --lst_start_hours and --lst_end_hours")
            sys.exit(-1)
#         uxtime = options.gps + 315964783
      
      if uxtime is not None :
         ( ra , dec , options.azim_deg , alt , options.za_deg ) = fits_beam.radec2azh( options.pointing_ra_deg, options.pointing_dec_deg, uxtime )
         print("INFO : converted (RA,DEC) = (%.4f,%.4f) [deg] into (AZ,ZA) = (%.4f,%.4f) [deg] for uxtime = %.2f" % (options.pointing_ra_deg, options.pointing_dec_deg, options.azim_deg, options.za_deg, uxtime ))
   
   print("###############################################################")
   print("PARAMATERS : ")
   print("###############################################################")
   print("Do plotting                = %s" % (options.do_plot))
   if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
      print("Pointing direction (ra,dec) = (%.4f,%.4f) [deg]" % (options.pointing_ra_deg,options.pointing_dec_deg))
   else :
      print("Pointing direction in RADEC not specified")   
      
   if options.azim_deg is not None and options.za_deg is not None : 
      print("Pointing direction (az,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg))
   else :
      print("Pointing direction in AZEL not specified")
   if options.lst_hours is not None :     
      print("Specified LST time         = %.4f [deg]" % (options.lst_hours))
   if options.freq_mhz is not None :
      print("Frequency                  = %.4f [MHz]" % (options.freq_mhz))
   else :
      print("Frequency                  = None") 
   if options.lst_start_hours is not None and options.lst_end_hours is not None :      
      print("LST range                  = %.4f - %.4f [MHz]" % (options.lst_start_hours,options.lst_end_hours))
   else :
      print("LST range not specified")
   if options.unixtime_start is not None and options.unixtime_end is not None :
      print("Unix time range            = %.2f - %.2f ( interval = %.2f seconds )" % (options.unixtime_start,options.unixtime_end,options.unixtime_interval))
   else :
      print("Unix time range not specified")
   if options.ut_start is not None and options.ut_end is not None :
      print("UTC time range = %s - %s ( interval = %.2f seconds )" % (options.ut_start,options.ut_end,options.unixtime_interval))
   else :
      print("UTC time range not specified")
   print("Debug level = %d" % (debug_level))
   print("###############################################################")

   return (options, args)

def save_output_file( freq, aot, pol, out_file_base ) :
   outfile_pol = ( "%s_%s.txt" % (out_file_base,pol) ) 
   out_f = open( outfile_pol , "w" )
   out_f.write( code_version_line )
   
   header_line = "# Frequency[MHz]  A/T[m^2/K]\n"
   out_f.write( header_line )
   
   len = freq.shape[0]
   n_lines = 0
   for i in range(0,len) :
      line = "%.4f %.8f\n" % (freq[i],aot[i])
      
      out_f.write( line )
      n_lines += 1
   
   out_f.close()
   
   return n_lines

def save_sens_vs_lst_file( lst_x, aot_x, sefd_x, lst_y, aot_y, sefd_y, out_file_base ) :
   outfile = ( "%s.txt" % (out_file_base) ) 
   out_f = open( outfile , "w" )
   out_f.write( code_version_line )
   
   header_line = "#  LST[h] A/T_x[m^2/K] SEFD_x[Jy] A/T_y[m^2/K] SEFD_y[Jy]\n"
   out_f.write( header_line )
   
   len = lst_x.shape[0]
   n_lines = 0
   for i in range(0,len) :
      line = "%.4f %.6f %.6f %.6f %.6f\n" % (lst_x[i], aot_x[i], sefd_x[i], aot_y[i], sefd_y[i])
      
      out_f.write( line )
      n_lines += 1
   
   out_f.close()
   print("DEBUG : saved sensitivity vs. LST to output file %s" % (outfile))
   
   return (outfile)

def save_sens_vs_freq_file( freq_x, aot_x, sefd_x, freq_y, aot_y, sefd_y, out_file_base ) :
   outfile = ( "%s.txt" % (out_file_base) ) 
   out_f = open( outfile , "w" )
   out_f.write( code_version_line )
   
   header_line = "#  FREQ[MHz] A/T_x[m^2/K] SEFD_x[Jy] A/T_y[m^2/K] SEFD_y[Jy] A/T_i[m^2/K] SEFD_i[Jy]\n"
   out_f.write( header_line )
   
   len = freq_x.shape[0]
   n_lines = 0
   for i in range(0,len) :
      sefd_i = 0.5*math.sqrt( sefd_x[i]*sefd_x[i] + sefd_y[i]*sefd_y[i] )
      aot_i = (2*1380.00)/sefd_i
   
      line = "%.4f %.6f %.6f %.6f %.6f %.6f %.6f\n" % (freq_x[i], aot_x[i], sefd_x[i], aot_y[i], sefd_y[i], sefd_i, aot_i)
      
      out_f.write( line )
      n_lines += 1
   
   out_f.close()
   print("DEBUG : saved sensitivity vs. FREQ to output file %s" % (outfile))
   
   return (outfile)


 
if __name__ == "__main__":
    (options, args) = parse_options(1)
        
    if ((options.azim_deg is not None and options.za_deg is not None) or (options.pointing_ra_deg is not None and options.pointing_dec_deg is not None)) and options.lst_hours is not None :
       # plotting A/T vs. frequency for a given pointing direction and LST :
       if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
          print("(RA,DEC) = (%.4f,%.4f) [deg] and LST = %.4f [h] specified -> getting A/T vs. frequency data and creating A/T vs. frequency plot ..." % (options.pointing_ra_deg,options.pointing_dec_deg,options.lst_hours))
       else : 
          print("(Azim,Za) and LST specified -> getting A/T vs. frequency data and creating A/T vs. frequency plot ...")

       (out_freq_x,out_aot_x,out_sefd_x,out_freq_y,out_aot_y,out_sefd_y,out_freq_i,out_aot_i,out_sefd_i) = get_sensitivity_azzalst( options.azim_deg, options.za_deg, options.lst_hours, receiver_temp_file=options.receiver_temp_file, station=options.station_name, ra_deg=options.pointing_ra_deg, dec_deg=options.pointing_dec_deg )
       
       if (out_freq_x is not None and out_aot_x is not None) or (out_freq_y is not None and out_aot_y is not None ) :
          if options.output_file is not None :
             if out_freq_x is not None and out_aot_x is not None :
                save_output_file( out_freq_x,out_aot_x, "X" , options.output_file )
             
             if out_freq_y is not None and out_aot_y is not None :
                save_output_file( out_freq_y,out_aot_y, "Y" , options.output_file )

             if out_freq_i is not None and out_aot_i is not None :
                save_output_file( out_freq_i , out_aot_i, "I" , options.output_file )

          if options.do_plot :
             if options.azim_deg is not None and options.za_deg is not None :
                info = "LST = %.1f h , (azim,za) = (%.2f,%.2f) [deg]" % (options.lst_hours,options.azim_deg,options.za_deg)
             else :
                info = "LST = %.1f h , (ra,dec) = (%.2f,%.2f) [deg]" % (options.lst_hours,options.pointing_ra_deg,options.pointing_dec_deg)
                if options.object_name is not None :
                   info += "\nObject = %s" % (options.object_name)
                
             plot_sensitivity( out_freq_x,out_aot_x, out_freq_y, out_aot_y, output_file_base=options.output_file, freq_i=out_freq_i, aot_i=out_aot_i, info=info )
             
          # do plot AoT vs. Freq :
    elif options.lst_hours is not None and options.freq_mhz is not None :
       # plotting sensitivity map of the whole sky for a given frequnecy and pointing direction :
       print("LST and frequency specified -> creating sensitivity (A/T) map over the whole hemisphere") 
       
       out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (options.lst_hours,options.freq_mhz)
       (azim_x,za_x,aot_x,sefd_x, azim_y,za_y,aot_y,sefd_y, azim_i,za_i,aot_i,sefd_i, out_txt_filename_X, out_txt_filename_Y, out_txt_filename_I ) = get_sensitivity_map( options.freq_mhz, options.lst_hours, output_file_base=out_fitsname_base,
                                                                                                                                                     out_fitsname_base=out_fitsname_base,
                                                                                                                                                     receiver_temperature=options.receiver_temperature, station=options.station_name )
       
       if ( azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None ) or ( azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None ) :
#           if options.output_file is not None :
#              if azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None :
                  # save map of sensitivity to file 
              
#              if azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None :
                  # save map of sensitivity to file
                  
           # test :
           if options.do_plot :
              # out_fitsname_base
              
              # Ensuring same color bar scale :
              vmax_value = max( max(aot_x) , max(aot_y) , max(aot_i) )
              vmax_value = 1.0
              print("DEBUG : common vmax for X,Y,I maps = %.4f" % (vmax_value))
              
              out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (options.lst_hours,options.freq_mhz)
              plot_sensitivity_map( azim_x, za_x, aot_x , out_fitsname_base=out_fitsname_base, save_text_file = options.save_text_file, do_plot=True, freq_mhz=options.freq_mhz, lst_h=options.lst_hours, pol="X", vmax_value=vmax_value )
              
              out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (options.lst_hours,options.freq_mhz)
              plot_sensitivity_map( azim_y, za_y, aot_y , out_fitsname_base=out_fitsname_base , save_text_file = options.save_text_file, do_plot=True, freq_mhz=options.freq_mhz, lst_h=options.lst_hours, pol="Y", vmax_value=vmax_value )


              print("DEBUG : %d / %d / %d" % (len(aot_x),len(aot_y),len(aot_i)))             
              out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (options.lst_hours,options.freq_mhz)
              plot_sensitivity_map( azim_i, za_i, aot_i, out_fitsname_base=out_fitsname_base , save_text_file = options.save_text_file, do_plot=True, freq_mhz=options.freq_mhz, lst_h=options.lst_hours, pol="I", s=5, vmax_value=vmax_value )
                            

    elif options.unixtime_start is not None and options.unixtime_end is not None :
        print("Plotting sensitivity for a specified time range unix time %.2f - %.2f" % (options.unixtime_start,options.unixtime_end))

        if options.freq_mhz is not None :
           if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
              print("\tPlotting for specific frequency = %.2f MHz and pointing direction (ra,dec) = (%.4f,%.4f) [deg]" % (options.freq_mhz,options.pointing_ra_deg,options.pointing_dec_deg))
           else :
              print("\tPlotting for specific frequency = %.2f MHz and pointing direction (azim,elev) = (%.4f,%.4f) [deg]" % (options.freq_mhz,options.azim_deg,options.za_deg))
           
           
           if ( options.azim_deg is not None and options.za_deg is not None ) or ( options.pointing_ra_deg is not None and options.pointing_dec_deg is not None ) :
              (uxtime_x,aot_x,sefd_x, uxtime_y,aot_y,sefd_y, uxtime_i,aot_i,sefd_i) = get_sensitivity_timerange( options.azim_deg, options.za_deg, options.freq_mhz, options.unixtime_start, options.unixtime_end, 
                                                                                          time_step=options.step_seconds, station=options.station_name, receiver_temperature=options.receiver_temperature,
                                                                                          ra_deg=options.pointing_ra_deg, dec_deg=options.pointing_dec_deg  )

              if options.do_plot :
                 info = ""
                 if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
                    info = "Freq. = %.3f [MHz]\n(ra,dec) = (%.4f,%.4f) [deg]" % ( options.freq_mhz,options.pointing_ra_deg,options.pointing_dec_deg)
                    if options.object_name is not None :
                       info += "\nObject = %s" % (options.object_name)
                 else :
                    info = "Freq. = %.3f [MHz]\n(azim,za) = (%.2f,%.2f) [deg]" % ( options.freq_mhz,options.azim_deg,options.za_deg)
              
                 out_fitsname_base = "sensitivity_uxtimerange%.2f-%.2f_az%.4fdeg_za%.4fdeg_freq%06.2fMHz" % (options.unixtime_start,options.unixtime_end,options.azim_deg,options.za_deg,options.freq_mhz)          
                 plot_sensitivity_vs_time( uxtime_x, aot_x, uxtime_y, aot_y, options.unixtime_start, options.unixtime_end, options.azim_deg,\
                                           options.za_deg, options.freq_mhz, output_file_base=out_fitsname_base, uxtime_i=uxtime_i, aot_i=aot_i, info=info,\
                                           save_text_file = options.save_text_file  )
              
              
           else :
              print("\tPlotting full sky maps in time - NOT YET IMPLEMENTED")  
        else :
           print("\tPlotting for range of frequencies (NOT YET IMPLEMENTED)")

    elif options.lst_start_hours is not None and options.lst_end_hours is not None :
       print("SENS_vs_LST : Plotting sensitivity for a specified LST %.2f - %.2f" % (options.lst_start_hours,options.lst_end_hours))
       
       if options.freq_mhz is not None :
           print("\tSENS_vs_LST : Plotting for specific frequency = %.2f MHz" % (options.freq_mhz))
           
           if ( options.azim_deg is not None and options.za_deg is not None ) or  ( options.pointing_ra_deg is not None and options.pointing_dec_deg is not None ) :
              
              if options.pointing_ra_deg is not None and options.pointing_dec_deg is not None :
                 print("\tSENS_vs_LST RADEC : Plotting for pointing direction (ra,dec) = (%.4f,%.4f) [deg]" % (options.pointing_ra_deg,options.pointing_dec_deg))
                 (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y, lst_i,aot_i,sefd_i ) = get_sensitivity_radec_lstrange( options.pointing_ra_deg, options.pointing_dec_deg, options.freq_mhz, options.lst_start_hours, options.lst_end_hours, 
                                                                                 time_step=options.step_seconds, station=options.station_name, receiver_temperature=options.receiver_temperature )

                 if options.do_plot :                 
                    out_fitsname_base = "sensitivity_lstrange%.2f-%.2f_ra%.4fdeg_dec%.4fdeg_freq%06.2fMHz_X" % (options.lst_start_hours, options.lst_end_hours,options.pointing_ra_deg,options.pointing_dec_deg,options.freq_mhz)          
                    plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, options.lst_start_hours, options.lst_end_hours, -1000, -1000, options.freq_mhz, output_file_base=out_fitsname_base, lst_i=lst_i,aot_i=aot_i, ra_deg=options.pointing_ra_deg, dec_deg=options.pointing_dec_deg )

              elif options.azim_deg is not None and options.za_deg is not None :
                 print("\tSENS_vs_LST : Plotting for pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg))
                 (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y, lst_i,aot_i,sefd_i ) = get_sensitivity_lstrange( options.azim_deg, options.za_deg, options.freq_mhz, options.lst_start_hours, options.lst_end_hours, 
                                                                                 time_step=options.step_seconds, station=options.station_name, receiver_temperature=options.receiver_temperature )

                 if options.do_plot :                 
                    out_fitsname_base = "sensitivity_lstrange%.2f-%.2f_az%.4fdeg_za%.4fdeg_freq%06.2fMHz_X" % (options.lst_start_hours, options.lst_end_hours,options.azim_deg,options.za_deg,options.freq_mhz)          
                    plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, options.lst_start_hours, options.lst_end_hours, options.azim_deg, options.za_deg, options.freq_mhz, output_file_base=out_fitsname_base, lst_i=lst_i,aot_i=aot_i )

              else :
                 print("\tSENS_vs_LST : Plotting for pointing direction (ra,dec) = (%.4f,%.4f) [deg]" % (options.pointing_ra_deg,options.pointing_dec_deg))
                 printf("TO BE IMPLEMENTED SOON !!!")
                    
              
              
           else :
              print("\tPlotting full sky maps in time - NOT YET IMPLEMENTED")  
       else :
           print("\tPlotting for range of frequencies (NOT YET IMPLEMENTED)")
       
       
    # TODO :
    # option to plot map of the sky in zenith projection showing sensitivity across the whole sky
    #    help : https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
    #           https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # 
    # option to plot sensitivity at a given pointing direction and frequency and over a specified time range 
    
        
    