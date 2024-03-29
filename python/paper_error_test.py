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
#     b/ A/T vs. lst at a particular pointing direction and frequency at zenith
#        python ./sensitivity_db.py --freq_mhz=154.88 --lst_start=0.00 --lst_end=24.00 --azim_deg=0 --za_deg=0 --do_plot      
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

debug_level = 0 


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
                             receiver_temp_file=None ) :
        
    global debug_level
    
    
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
    szSQL = "SELECT MIN(ABS(lst-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f" %  (lst_hours,lst_hours,db_lst_resolution, za_deg, db_ang_res_deg, az_deg, db_ang_res_deg )
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
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f" %  (lst_hours,(min_lst_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg )
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
       dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00-za_deg_db) , az_deg, (90.00 - za_deg) )
       
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
           sefd = numpy.NaN
           if aot > 0 :
              sefd = (2*1380.00)/aot
        
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) )
              
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
           print("\t\tDB Record ignored due to angular distance too large")
 

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_freq_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_freq_x, out_sefd_y, out_freq_y )

    return ( numpy.array(out_freq_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_freq_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_freq_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )


def unixtime2lst( unixtime ) :
    utc_str = datetime.utcfromtimestamp( unixtime ).strftime('%Y-%m-%d %H:%M:%S')
    t_utc = Time( utc_str, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))

    lst_hours = t_utc.sidereal_time('apparent').value
    
    return lst_hours


# get sensitivity vs. unix time for a specific polarisation :
def get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=5.00, 
                             receiver_temperature=None) :

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
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s'" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol )
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
       dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00-za_deg_db) , az_deg, (90.00 - za_deg) )
       
       if dist_value_deg < min_angular_distance_deg :
          min_angular_distance_deg = dist_value_deg
          closest_gridpoint_za_deg = za_deg_db
          closest_gridpoint_az_deg = azim_deg_db

    if min_angular_distance_deg > 1000 :
       print("ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg))
       return (None,None,None,None,None,None)
       
       
    print("Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg))


    used_db_record_id = []
        
    for unixtime in range( int(ux_start), int(ux_end)+1, time_step ) :
       lst_hours = unixtime2lst( unixtime )
       
#       print("DEBUG : ux=%.2f -> lst=%.2f" % (unixtime,lst_hours))
    
       cur = conn.cursor()       
       szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,a_eff,t_sys,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE lst>=%.4f AND lst<=%.4f AND frequency_mhz>=%.4f AND frequency_mhz<=%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s'" %  \
                 ( (lst_hours-db_lst_resolution),(lst_hours+db_lst_resolution),\
                   (freq_mhz-(db_freq_resolution_mhz+0.01)), (freq_mhz+(db_freq_resolution_mhz+0.01)),\
                   za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol\
                  )
       if debug_level >= 2 :
          print("DEBUG SQL get_sensitivity_timerange_single_pol : %s" % (szSQL))
       cur.execute( szSQL )
       rows = cur.fetchall()
       
       min_lst_distance = db_lst_resolution
       best_rec_id = -1
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
                 
       if debug_level >= 2 or min_lst_distance < db_lst_resolution :
          print("DEBUG : closest in time is min_lst_diff = %.4f [h] (id = %d)" % (min_lst_distance,best_rec_id))


       last_lst = -1000        
       for row in rows :
           if debug_level > 0 : 
              print(row)
        
           id = int( row[0] )
           
           if id in used_db_record_id :
              # skip already used records :
              continue
           
           azim_deg_db = float( row[1] )
           za_deg_db   = float( row[2] )
           freq_mhz    = float( row[3] )
           pol         = row[4]
           lst_db      = float( row[5] )
           unixtime_db = float( row[6] ) # when filled 
           gpstime     = float( row[7] )
           aot         = float( row[8] )
           sefd        = (2*1380.00)/aot
           a_eff       = float( row[9] )
           t_sys       = float( row[10] )
           t_rcv       = float( row[11] )
           t_ant       = float( row[12] )
           array_type  = int( row[13] )
           timestamp   = row[14]
           creator     = row[15]
           code_version = row[16]
           
           if receiver_temperature is not None and receiver_temperature >= 0.00 :
              # in case trcv is provided in the external text file (like a config file) :
              t_rcv = receiver_temperature
              aot = a_eff / ( t_ant + t_rcv )
              sefd        = (2*1380.00)/aot
           
           print("TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,min_angular_distance_deg,az_deg,za_deg))
                
           if best_rec_id < 0 or id == best_rec_id : 
              ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) ) 
        
              if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
                  if debug_level >= 4 :
                     is_in_list = False
                     if id in used_db_record_id :
                        is_in_list = True
                     print("DEBUG_FINAL id = %d : %.4f %.4f %.4f , list count = %d , is_in_list = %s" % (id,unixtime,aot,sefd,len(used_db_record_id),is_in_list))
                     
                  out_uxtime.append( unixtime )
                  out_aot.append( aot )
                  out_sefd.append( sefd ) 
                  
                  used_db_record_id.append( id ) 
              else :
                 print("\t\tDB Record ignored due to angular distance too large")
             
           if out_f is not None :
              line = "%.4f %.4f\n" % (unixtime,t_ant) 
              out_f.write( line )
           
           last_lst = lst_db
 
    if out_f is not None :
       out_f.close()

    return (out_uxtime, out_aot, out_sefd)

# get sensitivity vs. LST for a specific polarisation :
def get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None) :

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
    szSQL = "SELECT MIN(ABS(frequency_mhz-%.4f)) FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol, lst_start, lst_end )
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
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<=%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,(min_freq_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol, lst_start, lst_end )
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
       dist_value_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00-za_deg_db) , az_deg, (90.00 - za_deg) )
       
       if dist_value_deg < min_angular_distance_deg :
          min_angular_distance_deg = dist_value_deg
          closest_gridpoint_za_deg = za_deg_db
          closest_gridpoint_az_deg = azim_deg_db

    if min_angular_distance_deg > 1000 :
       print("ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg))
       return (None,None,None,None,None,None)
       
       
    print("Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg))


    
    cur = conn.cursor()       
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f ORDER BY LST ASC" %  (freq_mhz,(min_freq_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol, lst_start, lst_end )
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
                
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) ) 
        
        if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
            out_lst.append( lst_db )
            out_aot.append( aot )
            out_sefd.append( sefd ) 
        else :
            print("\t\tDB Record ignored due to angular distance too large")
 
    return (out_lst, out_aot, out_sefd)



# get sensitivity vs. time :
def get_sensitivity_timerange( az_deg , za_deg , freq_mhz, ux_start, ux_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=5.00, 
                             receiver_temperature=None ) :
        
    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_uxtime_x, out_aot_x, out_sefd_x) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    (out_uxtime_y, out_aot_y, out_sefd_y) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz,
                                                                                  receiver_temperature=receiver_temperature )

    # calculate SEFD_I - if possible (same sizes of arrays):
    ( out_uxtime_i , out_sefd_i , out_aot_i ) = calc_sefd_i( out_sefd_x, out_uxtime_x, out_sefd_y, out_uxtime_y )

    return ( numpy.array(out_uxtime_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_uxtime_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y),
             numpy.array(out_uxtime_i), numpy.array(out_aot_i) , numpy.array(out_sefd_i) )

def get_sensitivity_lstrange( az_deg , za_deg , freq_mhz, lst_start, lst_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, 
                             receiver_temperature=None) :

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
                             receiver_temperature=None, debug=True ) :
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
   
    if conn is None :    
       print("ERROR : could not connect to database %s" % (dbname_file))
       return (None,None,None,None,None,None)



    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE lst between %.4f and %.4f AND ABS(za_deg-%.4f)<%.4f AND ABS(azim_deg-%.4f)<%.4f AND ABS(frequency_mhz-%.4f)<%.4f ORDER BY lst ASC" %  (lst_start_h,lst_end_h,za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, freq_mhz, freq_resolution_mhz )
    if debug :
       print("SQL : %s" % (szSQL))
    cur.execute( szSQL )
    rows = cur.fetchall()
 
    out_lst_x = []
    out_aot_x  = []
    out_sefd_x = []
    
    out_lst_y = []
    out_aot_y  = []
    out_sefd_y = []
    
    if debug :
       print("   Fetched %d rows" % (len(rows)))
        
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
                              fig_size_x=20, fig_size_y=10, info=None ) :
                     
   legend_location = "upper center" # "upper right"

   # conversion of unix time to UTC :                               
   x_utc = []
   y_utc = []
   for ux in uxtime_x :
      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      x_utc.append( utc )
            
   for ux in uxtime_y :
      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      y_utc.append( utc )
   
   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   fig = plt.gcf()
   
   if uxtime_x is not None and aot_x is not None :
      ax_x =  plt.plot( x_utc, aot_x, point_x )
      # ax_x.tick_params(which='major', length=10, width=20, direction='inout')
#      ax_x.set_ylim(( min_ylimit,  max_ylimit ))
   
   if uxtime_y is not None and aot_y is not None :
      ax_y = plt.plot( y_utc, aot_y, point_y )

   if uxtime_i is not None and aot_i is not None :
      ax_i = plt.plot( x_utc, aot_i, point_i )
      
      max_i = max(aot_i)
      if max_i > max_ylimit :
         max_ylimit = max_i*1.1 # max_i + 10% max_i
#      ax_y.set_ylim(( min_ylimit,  max_ylimit ))

   # legend :
   if uxtime_x is not None and aot_x is not None and uxtime_y is not None and aot_y is not None :
      if uxtime_i is not None and aot_i is not None :
         plt.legend(('X polarisation','Y polarisation','Stokes I'), loc=legend_location, fontsize=20)
      else :
         plt.legend(('X polarisation','Y polarisation'), loc=legend_location,  fontsize=20 )
   else :
      if uxtime_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc=legend_location,  fontsize=20 )
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=legend_location,handles=[uxtime_y, uxtime_y])
   
   plt.xlabel('Time [UTC]' , fontsize=20 )
   plt.ylabel('Sensitivity A / T [m$^2$/K]' , fontsize=20 )
   plt.gcf().autofmt_xdate()
   
   ax_list = fig.axes
   ticks = plt.xticks
   x_start, x_end = ax_list[0].get_xlim()
   print("DEBUG : x_start = %.8f , x_end = %.8f" % (x_start,x_end))
#   ax_list[0].xaxis.set_ticks(numpy.arange(x_start, x_end, 10.00))
#   ax_list[0].xaxis.set_major_locator(md.MinuteLocator(interval=1))
#   ax_list[0].xaxis.set_major_locator(md.MinuteLocator(interval=10))
   if ( x_end - x_start ) >= 20.00 :
      ax_list[0].xaxis.set_ticks(numpy.arange( x_start, x_end, 4))
   else :
      ax_list[0].xaxis.set_ticks(numpy.arange( x_start, x_end, 2))
   
   
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M:%S')
#   ax.xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))
   max_ylimit0 = max_ylimit
   max_ylimit = max_ylimit*1.1 # + 10%
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
   if info is not None :
      # place a text box in upper left in axes coords
      min_x = min( x_utc )
      text( x_start, max_ylimit0, info , fontsize=20 ) # , transform=ax.transAxes, verticalalignment='top', bbox=props)
   
   plt.grid()
   
   if output_file_base is not None :
      outfile = ( "%s.png" % (output_file_base) )
      plt.savefig( outfile )
   
   plt.show()
   
def plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y,  lst_start, lst_end, azim_deg, za_deg, freq_mhz,
                              output_file_base=None, point_x='go', point_y='rx',
                              min_ylimit=0.00, max_ylimit=2.00,
                              do_show = True, save_output_path="./",
                              lst_i=None, aot_i=None, point_i='b+' ) :

   global web_interface_initialised
   
   plt.figure()
   if lst_x is not None and aot_x is not None :
      ax_x =  plt.plot( lst_x, aot_x, point_x )
   
   if lst_y is not None and aot_y is not None :
      ax_y = plt.plot( lst_y, aot_y, point_y )

   if lst_i is not None and aot_i is not None :
      ax_i = plt.plot( lst_i, aot_i, point_i )
      
      max_i = max( aot_i )
      if max_i > max_ylimit :
         max_ylimit = max_i * 1.1 # max_i + 10% of max_i

   # legend :
   if lst_x is not None and aot_x is not None and lst_y is not None and aot_y is not None :
      if lst_i is None and aot_i is None :      
         plt.legend(('X polarisation','Y polarisation'), loc='upper right')
      else :
         plt.legend(('X polarisation','Y polarisation','Stokes I'), loc='upper right')
   else :
      if lst_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc='upper right')
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[lst_y, lst_y])
   
   plt.xlabel('Local sidereal time [hours]')
   plt.ylabel('Sensitivity A / T [m$^2$/K]')
   plt.gcf().autofmt_xdate()
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M:%S')
#   ax.xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
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
                      fig_size_x=20, fig_size_y=10 ):

   global web_interface_initialised
   
   max_aot = 0.00
   min_freq = min(freq_x)
   
   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   ax = None
   if freq_x is not None and aot_x is not None :
      plt.plot( freq_x, aot_x, point_x )
      
      if max(aot_x) > max_aot :
         max_aot = max(aot_x)
   
   if freq_y is not None and aot_y is not None :
      plt.plot( freq_y, aot_y, point_y )
      
      if max(aot_y) > max_aot :
         max_aot = max(aot_y)


   if freq_i is not None and aot_i is not None :
      plt.plot( freq_i, aot_i, point_i )

      if max(aot_i) > max_aot :
         max_aot = max(aot_i)


   # legend :
   if freq_x is not None and aot_x is not None and freq_y is not None and aot_y is not None :
      if freq_i is not None and aot_i is not None :
         plt.legend(('X polarisation','Y polarisation','Stokes I'), loc='upper right' , fontsize=20)    
      else :
         plt.legend(('X polarisation','Y polarisation'), loc='upper right' , fontsize=20)
   else :
      if freq_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc='upper right' , fontsize=20)
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[freq_y, aot_y],fontsize=20)
   
   plt.xlabel('Frequency (MHz)' , fontsize=20 )
   plt.ylabel('Sensitivity A/T [m$^2$/K]' , fontsize=20 )
   plt.grid()

#   if max_aot < 0.5 :
#      max_ylimit = 0.5
      
   max_ylimit0 = max_ylimit
   max_ylimit = max_ylimit*1.1 # + 10%
   
   plt.ylim(( min_ylimit,  max_ylimit )) 
   
   if info is not None :
      # place a text box in upper left in axes coords
      text( min_freq, max_ylimit0*1.05, info , fontsize=20 ) # , transform=ax.transAxes, verticalalignment='top', bbox=props)

   
   outfile = None
   if output_file_base is not None :
      outfile = ( "%s.png" % (output_file_base) )
      plt.savefig( outfile )

      print("Saved output image to file %s" % (outfile))
      
      if web_interface_initialised :
         buf = BytesIO()
         plt.savefig( buf, format="png")
         # WORKS : see https://stackoverflow.com/questions/52368870/display-matplotlib-image-on-html-page-using-django         
         return (outfile,buf)
   
   plt.show()
   
   return( outfile , None )

def plot_sensitivity_map( azim_deg, za_deg, aot, out_fitsname_base="sensitivity" , save_text_file=False, do_plot=False, freq_mhz=0.00, lst_h=0.00, pol="Unknown", out_dir="./", s=3.5 ) :
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
      figsize=8
      fig=plt.figure(figsize=(figsize,0.6*figsize),dpi=300)
      axis('on')
      ax1=fig.add_subplot(1,1,1,polar=False)
    
      axis('off')
      ax2=fig.add_subplot(1,1,1,polar=True, frameon=False)
      ax2.set_theta_zero_location("N")
      ax2.set_theta_direction(-1)
      ax2.patch.set_alpha(0.0)
      ax2.tick_params(color='0.5', labelcolor='0.5')
      for spine in ax2.spines.values():
          spine.set_edgecolor('0.5')
      ax2.grid(which='major', color='0.5') 
    
      #Beamsky example:
#      vmax = 1.00
      vmax = max_sens
      im=ax1.imshow( sensitivity , interpolation='none', vmax=vmax )

      #Add colorbar on own axis
      cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
      cbar=fig.colorbar(im, cax=cbar_ax)
      cbar.set_label("A / T [m$^2$/K]")
 
      title = "Sensitivity map at frequency = %.2f MHz at lst = %.2f [h] (%s)" % (freq_mhz,lst_h,pol)
      ax1.set_title( title , fontsize=8 )
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
   parser.add_option('-l','--lst','--lst_hours',dest="lst_hours",default=None, help="Local sidreal time in hours [default %default]",metavar="float",type="float")

   # specific frequency in MHz 
   parser.add_option('-f','--freq_mhz','--frequency_mhz',dest="freq_mhz",default=None, help="Specific frequency in MHz [default %default]",metavar="float",type="float")

   # specify lst range :
   parser.add_option('--lst_start','--lst_start_hours',dest="lst_start_hours",default=None, help="Start time in Local sidreal time in hours [default %default]",metavar="float",type="float")
   parser.add_option('--lst_end','--lst_end_hours',dest="lst_end_hours",default=None, help="End time in Local sidreal time in hours [default %default]",metavar="float",type="float")

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
   
   # paper test :
   parser.add_option('--paper_test',dest="paper_test",default="freq", help="Type of paper test [default %default] other options lst, azza")

   (options, args) = parser.parse_args(sys.argv[idx:])
   
   
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
   
   print("###############################################################")
   print("PARAMATERS : ")
   print("###############################################################")
   print("paper_test = %s" % (options.paper_test))
   print("###############################################################")

   return (options, args)

def save_output_file( freq, aot, pol, out_file_base ) :
   outfile_pol = ( "%s_%s.txt" % (out_file_base,pol) ) 
   out_f = open( outfile_pol , "w" )
   
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
        
#    def get_sensitivity_azzalst( az_deg , za_deg , lst_hours , 
#                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
#                             db_lst_resolution=0.5, db_ang_res_deg=5.00,
#                             receiver_temp_file=None ) :

    uxtime = [1578443400,1578445200,1578447000,1578448800,1578450600,1578452400,1578454200,1578456000,1578457800,1578459600,1578461400,1578463200,1578465000,1578466800,1578468600,1578470400,1578472200,1578474000,1578475800,1578477600,1578479400,1578481200,1578483000,1578484800,1578486600,1578488400,1578490200,1578492000,1578493800,1578495600,1578497400,1578499200,1578501000,1578502800,1578504600,1578506400,1578508200,1578510000,1578511800,1578513600,1578515400,1578517200,1578519000,1578520800,1578522600,1578524400,1578526200,1578528000]
    out_progress_f = open("progress.txt","w")

    if options.paper_test == "freq" :    
       out_f = open("diff_xy_freq.txt","w")
       
       out_f.write( "# REL_DIFF_AoT_X REL_DIFF_AoT_Y AoT_X[ch] AoT_X[ch+1] AoT_Y[ch] AoT_Y[ch+1] FREQ[ch] LST[h] AZIM[deg] ZA[deg]\n" )
       for lst_idx in range(0,48) :
          out_progress_f.write("lst index = %d\n" % (lst_idx))
        
          lst = lst_idx/2.00        
          for za in range(60,0,-5) :
             for az in range(0,360,5) :
                ( out_freq_x,out_aot_x,out_sefd_x,out_freq_y,out_aot_y,out_sefd_y,out_freq_i,out_aot_i,out_sefd_i ) = get_sensitivity_azzalst( az,za,lst,"AAVS2" )
                
                l = len(out_freq_x) 
                for i in range(0,l-1) :
                   diff_x = (out_aot_x[i] - out_aot_x[i+1])/((out_aot_x[i] + out_aot_x[i+1])*0.5)
                   diff_y = (out_aot_y[i] - out_aot_y[i+1])/((out_aot_y[i] + out_aot_y[i+1])*0.5)
                   print("%.8f %.8f" % (diff_x,diff_y))
                   line = ("%.8f %.8f %.8f %.8f %.8f %.8f %.4f %.4f %.4f %.4f\n" % (diff_x,diff_y,out_aot_x[i],out_aot_x[i+1],out_aot_y[i],out_aot_y[i+1],out_freq_x[i],lst,az,za))
                   out_f.write( line )
                   
       out_f.close()                   
    elif options.paper_test == "lst" :
       # def get_sensitivity_azzalstrange( az_deg , za_deg , freq_mhz, lst_start_h, lst_end_h , 
       #                     station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
       #                     db_lst_resolution=0.5, db_ang_res_deg=5.00, freq_resolution_mhz=5.00,
       #                     receiver_temperature=None ) :
       out_f = open("diff_xy_lst.txt","w")
       out_f.write( "# REL_DIFF_AoT_X REL_DIFF_AoT_Y AoT_X[lst] AoT_X[lst+1] AoT_Y[lst] AoT_Y[lst+1] FREQ[MHz] LST[i] AZIM[deg] ZA[deg]\n")
       for freq_mhz in range (50,360,10) :
          out_progress_f.write("freq = %.2f MHz" % (freq_mhz))
       
          for za in range(60,0,-5) :
             for az in range(0,360,5) :
                ( out_lst_x , out_aot_x , out_sefd_x, out_lst_y , out_aot_y , out_sefd_y ) = get_sensitivity_azzalstrange( az, za, freq_mhz, 0, 24.00, "AAVS2" )
                
                l = len(out_lst_x) 
                for i in range(0,l-1) :
                   diff_x = (out_aot_x[i] - out_aot_x[i+1])/((out_aot_x[i] + out_aot_x[i+1])*0.5)
                   diff_y = (out_aot_y[i] - out_aot_y[i+1])/((out_aot_y[i] + out_aot_y[i+1])*0.5)
                   print("%.8f %.8f" % (diff_x,diff_y))
                   line = ("%.8f %.8f %.8f %.8f %.8f %.8f %.4f %.4f %.4f %.4f\n" % (diff_x,diff_y,out_aot_x[i],out_aot_x[i+1],out_aot_y[i],out_aot_y[i+1],freq_mhz,out_lst_x[i],az,za))
                   out_f.write( line )
                   
       out_f.close()
    elif options.paper_test == "azza" :
# get sensitivity map for a given LST and freq_mhz
#def get_sensitivity_map( freq_mhz, lst_hours, 
#                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
#                             db_lst_resolution=0.5, db_freq_resolution_mhz=10.00, output_file_base=None,
#                             out_fitsname_base="sensitivity_map_", output_dir="./", receiver_temperature=None ) :

       # just to calculate indexes of next za value for the same azimuth 
       out_f2 = open("diff_azza_all.txt","w")
       # diff_x_za,
       # diff_y_za,
       # out_aot_x[i],
       # out_aot_y[i],
       # out_azim_x[i],
       # out_za_x[i],
       # out_aot_x[ next_idx ],
       # out_aot_y[ next_idx ],
       # freq_mhz,
       # lst
       out_f2.write("# DIFF_TYPE REL_DIFF_ZA_X REL_DIFF_ZA_Y AoT_X[azza] AoT_Y[azza] AZIM[i] ZA[i] AoT_X[next za/azim] AoT_Y[next za/azim] FREQ[MHz] LST[h]\n")
       
       ( out_azim_x, out_za_x, out_aot_x, out_sefd_x, out_azim_y, out_za_y, out_aot_y, out_sefd_y, out_azim_y, out_za_y, out_aot_i, out_sefd_i, out_txt_filename_X, out_txt_filename_Y, out_txt_filename_I ) = get_sensitivity_map( 150.00, 0.00, "AAVS2" )
       len_az = len(out_azim_x)
       len_za = len(out_za_x)
       
       if len_az != len_za :
          printf("ERROR : len_az != len_za ( %d != %d )" % (len_az,len_za))
          sys.exit(-1)

       print("DEBUG : dimensions (%d,%d)" % (len_az,len_za))
          
       next_za_index = numpy.ones( len_za )*(-1)
       for i in range(0,len_za ):
          az_val = out_azim_x[i]
          za_val = out_za_x[i]
          
          for k in range(i+1,len_za ):
             print("DEBUG : looking for (az,za) = (%.2f,%.2f) [deg] , current value (%.2f,%.2f) [deg]" % (az_val,(za_val+5.00),out_azim_x[k],out_za_x[k]))
             if math.fabs(out_azim_x[k]-az_val)<=0.0001 and math.fabs(out_za_x[k]-(za_val+5.00))<=0.0001 :
                print("   FOUND at k = %d -> next_za_index[%d] = %d!!!" % (k,i,k))
                next_za_index[i] = k
                break

          if next_za_index[i] < 0 :
             print("ERROR : could not find next ZA value for i = %d at (az,za) = (%.2f,%.2f) [deg]" % (i,az_val,za_val))

       for i in range(0,len_za ):
          if out_za_x[i] < 84 : 
             if next_za_index[i] < 0 :
                print("ERROR : no next za index at i = %d at (az,za) = (%.2f,%.2f) [deg]-> exiting now" % (i,out_azim_x[i],out_za_x[i]))
                sys.exit(-1)


#       out_f = open("diff_xy_azza.txt","w")
       for freq_mhz in range (50,360,10) :
          for lst_idx in range(0,48) :       
             out_progress_f.write("freq = %.2f MHz , lst_idx = %d" % (freq_mhz,lst_idx))
             lst = lst_idx/2.00             

             outfile = ( "diff_xy_azza_%dMHz_lst%.1fh.txt" % (freq_mhz,lst) )
             out_f = open( outfile , "w" )
             
             ( out_azim_x, out_za_x, out_aot_x, out_sefd_x, out_azim_y, out_za_y, out_aot_y, out_sefd_y, out_azim_y, out_za_y, out_aot_i, out_sefd_i, out_txt_filename_X, out_txt_filename_Y, out_txt_filename_I ) = get_sensitivity_map( freq_mhz, lst, "AAVS2" )
             
             
             for i in range(0,len_az-1) :
                if out_za_x[i] <= 61 :             
                   diff_x_az = out_aot_x[i] - out_aot_x[i+1]
                   diff_y_az = out_aot_y[i] - out_aot_y[i+1]
                   
                   line = ("%.8f %.8f AZIM_DIFF\n" % (diff_x_az,diff_y_az))
                   out_f.write( line )
                   
                   if next_za_index[i] >= 0 :
                      next_idx = int(next_za_index[i])
                      # print("DEBUG : index = %d, value = %.2f" % (next_idx,out_aot_x[next_idx]))
                      diff_x_za = ( out_aot_x[i] - out_aot_x[ next_idx ] ) / (0.5*(out_aot_x[i] + out_aot_x[ next_idx ]))
                      diff_y_za = ( out_aot_y[i] - out_aot_y[ next_idx ] ) / (0.5*(out_aot_y[i] + out_aot_y[ next_idx ]))
                   
                      line = ("%.8f %.8f ZA_DIFF %.8f %.8f at ( %.8f , %.8f ) [deg] | %.8f %.8f at ( %.8f ,%.8f ) [deg] | %.8f %.8f\n" % (diff_x_za,diff_y_za,out_aot_x[i],out_aot_y[i],out_azim_x[i],out_za_x[i],out_aot_x[ next_idx ],out_aot_y[ next_idx ],out_azim_x[next_idx],out_za_x[next_idx],freq_mhz,lst))
                      out_f.write( line )
                      
                      line = ("ZA_DIFF %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n" % (diff_x_za,diff_y_za,out_aot_x[i],out_aot_y[i],out_azim_x[i],out_za_x[i],out_aot_x[ next_idx ],out_aot_y[ next_idx ],freq_mhz,lst))
                      out_f2.write( line )
                                                                                                     

       
             out_f.close()

       out_f2.close()      
    else :
       print("ERROR : not implemented yet !")
                
             
       
    out_progress_f.close()
    
    