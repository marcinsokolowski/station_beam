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

# script for quering SQLITE3 or PostgreSQL databases for sensitivity values :
# HELP : python : https://www.sqlitetutorial.net/sqlite-python/ , https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/
#        Plot vs. time : /home/msok/ska/aavs/aavs0.5/trunk/analysis/MWAdas/scripts$ plot_delay.py
# 
# Google : "python interpolation on sphere"
#        https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.interpolate.SmoothSphereBivariateSpline.html
#        https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.LSQSphereBivariateSpline.html
#        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.SmoothSphereBivariateSpline.html


def create_connection_sqlite3( db_file ):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
 
    return conn

def calc_anglular_distance_degrees( azim1_deg, za1_deg, azim2_deg , za2_deg ) :
    azim1_rad = azim1_deg * (math.pi/180.00)
    elev1_rad = ( 90.00 - za1_deg ) * (math.pi/180.00)

    azim2_rad = azim2_deg * (math.pi/180.00)
    elev2_rad = ( 90.00 - za2_deg ) * (math.pi/180.00)

    cos_value = math.sin(elev1_rad)*math.sin(elev2_rad) + math.cos(elev1_rad)*math.cos(elev2_rad)*math.cos( azim1_rad - azim2_rad )
    dist_value_deg = math.acos( cos_value )*(180.00/math.pi)
    
    return dist_value_deg


# get sensitivity vs. frequency for given pointing direction [degrees] and lst [hours]:
def get_sensitivity_azzalst( az_deg , za_deg , lst_hours , 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00 ) :
        
    global debug_level
    
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
  
    # select closest LST :
    cur = conn.cursor()
    szSQL = "SELECT MIN(ABS(lst-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f" %  (lst_hours,lst_hours,db_lst_resolution, za_deg, db_ang_res_deg, az_deg, db_ang_res_deg )
    print "DEBUG SQL1 : %s" % (szSQL)
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_lst_distance = None
    for row in rows:
       min_lst_distance = float( row[0] )
    print "Best LST in database is closer than %.8f [hours]" % (min_lst_distance)
    
    if min_lst_distance is None :
       print "ERROR no records in the database exist closer than %.4f hours in LST in the direction of (azim,za) = (%.4f,%.4f) [deg]" % (db_lst_resolution,az_deg,za_deg)
       return (None,None,None,None,None,None)
 
    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f" %  (lst_hours,(min_lst_distance+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg )
    print "DEBUG SQL2 : %s" % (szSQL)
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
       print "ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg)
       return (None,None,None,None,None,None)
       
       
    print "Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg)
     
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
        sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) )
              
        print "TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,ang_distance_deg,az_deg,za_deg)
        
        
        
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
               print "ERROR : unknown polarisation = %s" % (pol)
        else :
           print "\t\tDB Record ignored due to angular distance too large"
 

    return ( numpy.array(out_freq_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_freq_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )


def unixtime2lst( unixtime ) :
    utc_str = datetime.utcfromtimestamp( unixtime ).strftime('%Y-%m-%d %H:%M:%S')
    t_utc = Time( utc_str, scale='utc', location=(MWA_POS.lon.value, MWA_POS.lat.value ))

    lst_hours = t_utc.sidereal_time('apparent').value
    
    return lst_hours


# get sensitivity vs. unix time for a specific polarisation :
def get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, ) :

    out_uxtime = []
    out_aot    = []
    out_sefd   = []

    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
  
    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s'" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol )
    print "DEBUG SQL2 : %s" % (szSQL)
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
       print "ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg)
       return (None,None,None,None,None,None)
       
       
    print "Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg)


    
    for unixtime in range( int(ux_start), int(ux_end)+1, time_step ) :
       lst_hours = unixtime2lst( unixtime )
    
       cur = conn.cursor()       
       szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND ABS(lst-%.4f)<%.4f and polarisation='%s'" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, lst_hours, db_lst_resolution, pol )
       if debug_level >= 2 :
          print "DEBUG SQL get_sensitivity_timerange_single_pol : %s" % (szSQL)
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
          print "DEBUG : closest in time is min_lst_diff = %.4f [h] (id = %d)" % (min_lst_distance,best_rec_id)

        
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
           sefd        = (2*1380.00)/aot
           a_eff       = float( row[9] )
           t_rcv       = float( row[10] )
           t_ant       = float( row[11] )
           array_type  = int( row[12] )
           timestamp   = row[13]
           creator     = row[14]
           code_version = row[15]
           
           print "TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,min_angular_distance_deg,az_deg,za_deg)
                
           if best_rec_id < 0 or id == best_rec_id : 
              ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) ) 
        
              if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
                  out_uxtime.append( unixtime )
                  out_aot.append( aot )
                  out_sefd.append( sefd ) 
              else :
                 print "\t\tDB Record ignored due to angular distance too large"
 
    return (out_uxtime, out_aot, out_sefd)

# get sensitivity vs. LST for a specific polarisation :
def get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, ) :

    out_lst = []
    out_aot    = []
    out_sefd   = []

    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
  
    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol, lst_start, lst_end )
    print "DEBUG SQL2 : %s" % (szSQL)
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
       print "ERROR : no record in DB close enough to the requested pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (az_deg,za_deg)
       return (None,None,None,None,None,None)
       
       
    print "Closest pointing direction in DB is at (azim,za) = (%.4f,%.4f) [deg] in angular distance = %.4f [deg]" % (closest_gridpoint_az_deg,closest_gridpoint_za_deg,min_angular_distance_deg)


    
    cur = conn.cursor()       
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(frequency_mhz-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (azim_deg-%.4f)<%.4f AND polarisation='%s' AND lst>%.2f AND lst<%.2f ORDER BY LST ASC" %  (freq_mhz,(db_freq_resolution_mhz+0.01), za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, pol, lst_start, lst_end )
    if debug_level >= 2 :
       print "DEBUG SQL get_sensitivity_timerange_single_pol : %s" % (szSQL)
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
        sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
           
        print "TEST : %d , freq_mhz = %.4f [MHz] , (azim_deg_db,za_deg_db) = (%.4f,%.4f) [deg] in %.8f [deg] distance from requested (azim_deg,za_deg) = (%.4f,%.4f) [deg]" % (id,freq_mhz,azim_deg_db,za_deg_db,min_angular_distance_deg,az_deg,za_deg)
                
        ang_distance_deg = calc_anglular_distance_degrees( azim_deg_db, (90.00 - za_deg_db) , az_deg, (90.00 - za_deg) ) 
        
        if ang_distance_deg <= (min_angular_distance_deg+0.01) :        
            out_lst.append( lst_db )
            out_aot.append( aot )
            out_sefd.append( sefd ) 
        else :
            print "\t\tDB Record ignored due to angular distance too large"
 
    return (out_lst, out_aot, out_sefd)



# get sensitivity vs. time :
def get_sensitivity_timerange( az_deg , za_deg , freq_mhz, ux_start, ux_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, ) :
        
    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_uxtime_x, out_aot_x, out_sefd_x) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz )

    (out_uxtime_y, out_aot_y, out_sefd_y) = get_sensitivity_timerange_single_pol( az_deg , za_deg , freq_mhz, ux_start, ux_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz )

    return ( numpy.array(out_uxtime_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_uxtime_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )

def get_sensitivity_lstrange( az_deg , za_deg , freq_mhz, lst_start, lst_end, time_step=300,
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, db_freq_resolution_mhz=10.00, ) :

    global debug_level
                                 

# Change to loop over time range in steps of step      
    (out_lst_x, out_aot_x, out_sefd_x) = get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol='X', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz )

    (out_lst_y, out_aot_y, out_sefd_y) = get_sensitivity_lstrange_single_pol( az_deg , za_deg , freq_mhz, lst_start, lst_end, pol='Y', time_step=time_step,
                                                                                  station=station, db_base_name=db_base_name, db_path=db_path, 
                                                                                  db_lst_resolution=db_lst_resolution, db_ang_res_deg=db_ang_res_deg, db_freq_resolution_mhz=db_freq_resolution_mhz )

    return ( numpy.array(out_lst_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_lst_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )


# get sensitivity map for a given LST and freq_mhz
def get_sensitivity_map( freq_mhz, lst_hours, 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_freq_resolution_mhz=10.00, output_file_base=None ) :
                       
    global debug_level                        
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )
  
    # select closest LST :
    cur = conn.cursor()
    szSQL = "SELECT MIN(ABS(lst-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<=%.4f AND ABS(frequency_mhz-%.4f)<=%.4f" %  (lst_hours,lst_hours,db_lst_resolution, freq_mhz, db_freq_resolution_mhz )
    print "DEBUG SQL_best_lst : %s" % (szSQL)
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_lst_distance = None
    for row in rows:
       min_lst_distance = float( row[0] )
    print "Best LST in database is closer than %.8f [hours]" % (min_lst_distance)
    
    if min_lst_distance is None :
       print "ERROR no records in the database exist closer than %.4f hours in LST in the direction of (azim,za) = (%.4f,%.4f) [deg]" % (db_lst_resolution,az_deg,za_deg)
       return (None,None,None,None,None,None)

    # select closest FREQUENCY :
    cur = conn.cursor()
    szSQL = "SELECT MIN(ABS(frequency_mhz-%.8f)) FROM Sensitivity WHERE ABS(lst-%.4f)<=%.4f AND ABS(frequency_mhz-%.4f)<=%.4f" %  (freq_mhz, lst_hours, min_lst_distance, freq_mhz, db_freq_resolution_mhz )
    print "DEBUG SQL_best_freq : %s" % (szSQL)
    cur.execute( szSQL )
    rows = cur.fetchall()
    min_freq_distance = None
    for row in rows:
       print(row)
       min_freq_distance = float( row[0] )
    print "Best FREQ in database is closer than %.8f [MHz]" % (min_freq_distance)
    
    if min_freq_distance is None :
       print "ERROR no records in the database exist closer than %.4f MHz in FREQ at LST around %.4f [hours]" % (db_freq_resolution_mhz,lst_hours)
       return (None,None,None,None,None,None)
 
    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<=%.8f AND ABS(frequency_mhz-%.4f)<=%.8f ORDER BY za_deg, azim_deg ASC" %  (lst_hours,(min_lst_distance+0.01), freq_mhz, (min_freq_distance+0.1) )
    print "DEBUG SQL_main : %s" % (szSQL)
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

    if output_file_base :
       out_txt_filename_X = out_fitsname_base + "_X.txt"
       out_txt_filename_Y = out_fitsname_base + "_Y.txt"
       
       out_txt_x_f = open( out_txt_filename_X , "w" )
       out_txt_y_f = open( out_txt_filename_Y , "w" )
       
       line = "# AZIM[deg] ZA[deg] A/T[m^2/K]"
       out_txt_x_f.write( line + "    for X polarisation\n" )
       out_txt_y_f.write( line + "    for Y polarisation\n" )
        
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
        sefd        = (2*1380.00)/aot
        a_eff       = float( row[9] )
        t_rcv       = float( row[10] )
        t_ant       = float( row[11] )
        array_type  = int( row[12] )
        timestamp   = row[13]
        creator     = row[14]
        code_version = row[15]
        
        if debug_level >= 2 :
           print "TEST : %d , freq_mhz = %.4f [MHz]" % (id,freq_mhz)
        
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
           print "ERROR : unknown polarisation = %s" % (pol)
 
    if out_txt_x_f is not None :
       out_txt_x_f.close()

    if out_txt_y_f is not None :
       out_txt_y_f.close()
       

    return ( numpy.array(out_azim_x), numpy.array(out_za_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_azim_y), numpy.array(out_za_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )



# 
def get_sensitivity_azzalstrange( az_deg , za_deg , freq_mhz, lst_start_h, lst_end_h , 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00, freq_resolution_mhz=5.00 ) :
                             
    # connect to the database :                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)       
    conn = create_connection_sqlite3( dbname_file )


    # get requested data :
    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE lst between (%.4f and %.4f) AND (za_deg-%.4f)<%.4f AND (az_deg-%.4f)<%.4f AND ABS(frequency_mhz-%.4f)<%.4f ORDER BY lst ASC" %  (lst_start_h,lst_end_h,za_deg, db_ang_res_deg, az_deg, db_ang_res_deg, freq_mhz, freq_resolution_mhz )
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
        
        print "TEST : %d , freq_mhz = %.4f [MHz]" % (id,freq_mhz)
        
        if pol == "X" :
           out_lst_x.append( lst_db )
           out_aot_x.append( aot )
           out_sefd_x.append( sefd ) 
        elif pol == "Y" :
           out_lst_y.append( lst_db )
           out_aot_y.append( aot )
           out_sefd_y.append( sefd )            
        else :
           print "ERROR : unknown polarisation = %s" % (pol)
 

    return ( numpy.array(out_lst_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_lst_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )



def plot_sensitivity_vs_time( uxtime_x, aot_x, uxtime_y, aot_y,  unixtime_start, unixtime_end, azim_deg, za_deg, freq_mhz,
                              output_file_base=None, point_x='go', point_y='rx',
                              min_ylimit=0.00, max_ylimit=2.00 ) :
                     
   # conversion of unix time to UTC :                               
   x_utc = []
   y_utc = []
   for ux in uxtime_x :
      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      x_utc.append( utc )
   for ux in uxtime_y :
      utc = datetime.utcfromtimestamp(ux).strftime('%Y-%m-%dT%H:%M')
      y_utc.append( utc )
   
   plt.figure()
   if uxtime_x is not None and aot_x is not None :
      ax_x =  plt.plot( x_utc, aot_x, point_x )
#      ax_x.set_ylim(( min_ylimit,  max_ylimit ))
   
   if uxtime_y is not None and aot_y is not None :
      ax_y = plt.plot( y_utc, aot_y, point_y )
#      ax_y.set_ylim(( min_ylimit,  max_ylimit ))

   # legend :
   if uxtime_x is not None and aot_x is not None and uxtime_y is not None and aot_y is not None :
      plt.legend(('X polarisation','Y polarisation'), loc='upper right')
   else :
      if uxtime_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc='upper right')
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[uxtime_y, uxtime_y])
   
   plt.xlabel('Time [UTC]')
   plt.ylabel('Sensitivity A/T [m^2/K]')
   plt.gcf().autofmt_xdate()
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M:%S')
#   ax.xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
   plt.grid()
   
   if output_file_base is not None :
      outfile = ( "%s.png" % (output_file_base) )
      plt.savefig( outfile )
   
   plt.show()
   
def plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y,  lst_start, lst_end, azim_deg, za_deg, freq_mhz,
                              output_file_base=None, point_x='go', point_y='rx',
                              min_ylimit=0.00, max_ylimit=2.00 ) :
                     
   plt.figure()
   if lst_x is not None and aot_x is not None :
      ax_x =  plt.plot( lst_x, aot_x, point_x )
   
   if lst_y is not None and aot_y is not None :
      ax_y = plt.plot( lst_y, aot_y, point_y )

   # legend :
   if lst_x is not None and aot_x is not None and lst_y is not None and aot_y is not None :
      plt.legend(('X polarisation','Y polarisation'), loc='upper right')
   else :
      if lst_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc='upper right')
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[lst_y, lst_y])
   
   plt.xlabel('Local sidereal time [hours]')
   plt.ylabel('Sensitivity A/T [m^2/K]')
   plt.gcf().autofmt_xdate()
#   ax=plt.gca()
#   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
#   xfmt = md.DateFormatter('%H:%M:%S')
#   ax.xaxis.set_major_formatter(xfmt)
#   ax.set_ylim(( min_ylimit,  max_ylimit ))
   plt.ylim(( min_ylimit,  max_ylimit ))  
   
   plt.grid()
   
   if output_file_base is not None :
      outfile = ( "%s.png" % (output_file_base) )
      plt.savefig( outfile )
   
   plt.show()
   

def plot_sensitivity( freq_x, aot_x, freq_y, aot_y, output_file_base=None, point_x='go', point_y='rx' ):

   plt.figure()
   if freq_x is not None and aot_x is not None :
      plt.plot( freq_x, aot_x, point_x )
   
   if freq_y is not None and aot_y is not None :
      plt.plot( freq_y, aot_y, point_y )

   # legend :
   if freq_x is not None and aot_x is not None and freq_y is not None and aot_y is not None :
      plt.legend(('X polarisation','Y polarisation'), loc='upper right')
   else :
      if freq_x is not None and aot_x is not None :
         plt.legend(('X polarisation'), loc='upper right')
      else :
         plt.legend(bbox_to_anchor=(0.68, 0.82),loc=3,handles=[freq_y, aot_y])
   
   plt.xlabel('Frequency (MHz)')
   plt.ylabel('Sensitivity A/T [m^2/K]')
   plt.grid()
   
   if output_file_base is not None :
      outfile = ( "%s.png" % (output_file_base) )
      plt.savefig( outfile )
   
   plt.show()

def plot_sensitivity_map( azim_deg, za_deg, aot, out_fitsname_base="sensitivity" ) :
   from scipy.interpolate import SmoothSphereBivariateSpline    

   azim_rad = azim_deg*(numpy.pi/180.0)
   za_rad   = za_deg*(numpy.pi/180.0)
   lut = SmoothSphereBivariateSpline( za_rad, azim_rad, aot , s=3.5)
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
   for i in range(0,za_rad.shape[0]) :
      for j in range(0,za_rad.shape[1]) :
         a_rad = az_rad[ i , j ]
         z_rad = za_rad[ i , j ]
         
         aot_value = lut( z_rad, a_rad )
         if aot_value >= 0 :
            # ignore non-physical negative values :
            sensitivity[ i , j ] = aot_value
         

   fits_beam.save_fits( sensitivity , out_fitsname_base + ".fits" )
   fits_beam.save_fits( az_rad*(180.00/numpy.pi) , out_fitsname_base + "_azim_deg.fits" )
   fits_beam.save_fits( za_rad*(180.00/numpy.pi) , out_fitsname_base + "_za_deg.fits" )
         
   # save txt file :
   if False : # saving of full (interpolated) map as in FITS file is turned of as the files are too large :
      out_txt_filename = out_fitsname_base + ".txt"
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
   
   return (lut,sensitivity)
 

 
def parse_options(idx):
   usage="Usage: %prog [options]\n"
   usage += "Different plots / data can be obtainted:\n"
   usage += " 1/ to plot AoT vs. freq. for a given pointing direction (az,za) and lst time use options --azim_deg, --za_deg, --lst_hours\n"
   usage += " 2/ to plot AoT vs. time for a specified frequency use options --freq_mhz, --lst_start, --lst_end\n"
   usage += " 3/ to plot sensitivity map over the entire sky at a specified LST and frequency specify options : --lst_hours and --freq_mhz\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # General parameters :
   parser.add_option('-s','--station_name','--station',dest="station_name",default="EDA2", help="Station name [default %default]")

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-p','--plot','--do_plot',action="store_true",dest="do_plot",default=False, help="Plot")
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
   
   print "###############################################################"
   print "PARAMATERS : "
   print "###############################################################"
   print "Do plotting                = %s" % (options.do_plot)
   if options.azim_deg is not None and options.za_deg is not None : 
      print "Pointing direction (az,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg)
   else :
      print "Pointing direction not specified"
   if options.lst_hours is not None :     
      print "Specified LST time         = %.4f [deg]" % (options.lst_hours)
   if options.freq_mhz is not None :
      print "Frequency                  = %.4f [MHz]" % (options.freq_mhz)
   else :
      print "Frequency                  = None" 
   if options.lst_start_hours is not None and options.lst_end_hours is not None :      
      print "LST range                  = %.4f - %.4f [MHz]" % (options.lst_start_hours,options.lst_end_hours)
   else :
      print "LST range not specified"
   if options.unixtime_start is not None and options.unixtime_end is not None :
      print "Unix time range            = %.2f - %.2f ( interval = %.2f seconds )" % (options.unixtime_start,options.unixtime_end,options.unixtime_interval)
   else :
      print "Unix time range not specified"
   if options.ut_start is not None and options.ut_end is not None :
      print "UTC time range = %s - %s ( interval = %.2f seconds )" % (options.ut_start,options.ut_end,options.unixtime_interval)
   else :
      print "UTC time range not specified"
   print "###############################################################"

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
 
if __name__ == "__main__":
    (options, args) = parse_options(1)

    if options.azim_deg is not None and options.za_deg is not None and options.lst_hours is not None :
       # plotting A/T vs. frequency for a given pointing direction and LST :
       (out_freq_x,out_aot_x,out_sefd_x,out_freq_y,out_aot_y,out_sefd_y) = get_sensitivity_azzalst( options.azim_deg, options.za_deg, options.lst_hours )
       
       if (out_freq_x is not None and out_aot_x is not None) or (out_freq_y is not None and out_aot_y is not None ) :
          if options.output_file is not None :
             if out_freq_x is not None and out_aot_x is not None :
                save_output_file( out_freq_x,out_aot_x, "X" , options.output_file )
             
             if out_freq_y is not None and out_aot_y is not None :
                save_output_file( out_freq_y,out_aot_y, "Y" , options.output_file )

          if options.do_plot :
             plot_sensitivity( out_freq_x,out_aot_x, out_freq_y, out_aot_y, output_file_base=options.output_file )
          # do plot AoT vs. Freq :
    elif options.lst_hours is not None and options.freq_mhz is not None :
       # plotting sensitivity map of the whole sky for a given frequnecy and pointing direction :
       out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (options.lst_hours,options.freq_mhz)
       (azim_x,za_x,aot_x,sefd_x,azim_y,za_y,aot_y,sefd_y) = get_sensitivity_map( options.freq_mhz, options.lst_hours, output_file_base=out_fitsname_base )
       
       if ( azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None ) or ( azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None ) :
#           if options.output_file is not None :
#              if azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None :
                  # save map of sensitivity to file 
              
#              if azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None :
                  # save map of sensitivity to file
                  
           # test :
           if options.do_plot :
              # out_fitsname_base
              out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz_X" % (options.lst_hours,options.freq_mhz)
              plot_sensitivity_map( azim_x, za_x, aot_x , out_fitsname_base=out_fitsname_base)
              
              out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz_Y" % (options.lst_hours,options.freq_mhz)
              plot_sensitivity_map( azim_y, za_y, aot_y , out_fitsname_base=out_fitsname_base)

    elif options.unixtime_start is not None and options.unixtime_end is not None :
        print "Plotting sensitivity for a specified time range unix time %.2f - %.2f" % (options.unixtime_start,options.unixtime_end)

        if options.freq_mhz is not None :
           print "\tPlotting for specific frequency = %.2f MHz" % (options.freq_mhz)
           
           if options.azim_deg is not None and options.za_deg is not None :
              print "\tPlotting for pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg)
              (uxtime_x,aot_x,sefd_x, uxtime_y,aot_y,sefd_y) = get_sensitivity_timerange( options.azim_deg, options.za_deg, options.freq_mhz, options.unixtime_start, options.unixtime_end, time_step=options.step_seconds, station=options.station_name )

              if options.do_plot :
                 out_fitsname_base = "sensitivity_uxtimerange%.2f-%.2f_az%.4fdeg_za%.4fdeg_freq%06.2fMHz_X" % (options.unixtime_start,options.unixtime_end,options.azim_deg,options.za_deg,options.freq_mhz)          
                 plot_sensitivity_vs_time( uxtime_x, aot_x, uxtime_y, aot_y, options.unixtime_start, options.unixtime_end, options.azim_deg, options.za_deg, options.freq_mhz, output_file_base=out_fitsname_base )
              
              
           else :
              print "\tPlotting full sky maps in time - NOT YET IMPLEMENTED"  
        else :
           print "\tPlotting for range of frequencies (NOT YET IMPLEMENTED)"

    elif options.lst_start_hours is not None and options.lst_end_hours is not None :
       print "SENS_vs_LST : Plotting sensitivity for a specified LST %.2f - %.2f" % (options.lst_start_hours,options.lst_end_hours)
       
       if options.freq_mhz is not None :
           print "\tSENS_vs_LST : Plotting for specific frequency = %.2f MHz" % (options.freq_mhz)
           
           if options.azim_deg is not None and options.za_deg is not None :
              print "\tSENS_vs_LST : Plotting for pointing direction (azim,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg)
              (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y) = get_sensitivity_lstrange( options.azim_deg, options.za_deg, options.freq_mhz, options.lst_start_hours, options.lst_end_hours, time_step=options.step_seconds, station=options.station_name )

              if options.do_plot :
                 out_fitsname_base = "sensitivity_lstrange%.2f-%.2f_az%.4fdeg_za%.4fdeg_freq%06.2fMHz_X" % (options.lst_start_hours, options.lst_end_hours,options.azim_deg,options.za_deg,options.freq_mhz)          
                 plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, options.lst_start_hours, options.lst_end_hours, options.azim_deg, options.za_deg, options.freq_mhz, output_file_base=out_fitsname_base )
              
              
           else :
              print "\tPlotting full sky maps in time - NOT YET IMPLEMENTED"  
       else :
           print "\tPlotting for range of frequencies (NOT YET IMPLEMENTED)"
       
       
    # TODO :
    # option to plot map of the sky in zenith projection showing sensitivity across the whole sky
    #    help : https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
    #           https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # 
    # option to plot sensitivity at a given pointing direction and frequency and over a specified time range 
    
        
    