# SCRIPT for getting station sensitivity from the sensitivity database (SQLITE3 or PostgreSQL) :
# example commands :
#    python ./sensitivity_db.py --azim_deg=0 --za_deg=0 --lst=15.4 --out_file="Azim0_Za0_Lst15.4hours" --do_plot 
#    Verify results against files in templates/ Azim0_Za0_Lst15.4hours_X.txt and Azim0_Za0_Lst15.4hours_Y.txt
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
import matplotlib.pyplot as plt

# local packages :
import beam_tools

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

# get sensitivity vs. frequency for given pointing direction [degrees] and lst [hours]:
def get_sensitivity_azzalst( az_deg , za_deg , lst_hours , 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00 ) :
                             
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
 
    out_freq_x = []
    out_aot_x  = []
    out_sefd_x = []
    
    out_freq_y = []
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
           out_freq_x.append( freq_mhz )
           out_aot_x.append( aot )
           out_sefd_x.append( sefd ) 
        elif pol == "Y" :
           out_freq_y.append( freq_mhz )
           out_aot_y.append( aot )
           out_sefd_y.append( sefd )            
        else :
           print "ERROR : unknown polarisation = %s" % (pol)
 

    return ( numpy.array(out_freq_x), numpy.array(out_aot_x) , numpy.array(out_sefd_x),
             numpy.array(out_freq_y), numpy.array(out_aot_y) , numpy.array(out_sefd_y) )


# get sensitivity map for a given LST and freq_mhz
def get_sensitivity_map( freq_mhz, lst_hours, 
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="sql/", 
                             db_lst_resolution=0.5, db_freq_resolution_mhz=10.00 ) :
                       
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
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<=%.8f AND ABS(frequency_mhz-%.4f)<=%.8f" %  (lst_hours,(min_lst_distance+0.01), freq_mhz, (min_freq_distance+0.1) )
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
        
        if pol == "X" :
           out_azim_x.append( azim_deg_db )
           out_za_x.append( za_deg_db )
           out_aot_x.append( aot )
           out_sefd_x.append( sefd ) 
        elif pol == "Y" :
           out_azim_y.append( azim_deg_db )
           out_za_y.append( za_deg_db )
           out_aot_y.append( aot )
           out_sefd_y.append( sefd )            
        else :
           print "ERROR : unknown polarisation = %s" % (pol)
 

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

def plot_sensitivity_map( azim_deg, za_deg, aot, out_fitsname="sensitivity.fits" ) :
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
   (az_deg,za_deg,x,y,z,d) = beam_tools.makeAZZA(npix=256,projection='SIN',azim_from_north=True, return_all=True )    
   
   
   sensitivity = numpy.zeros( az_deg.shape )
   for i in range(0,za_deg.shape[0]) :
      for j in range(0,za_deg.shape[1]) :
         az_rad = az_deg[ i , j ] * (numpy.pi/180.00)
         za_rad = za_deg[ i , j ] * (numpy.pi/180.00)
         
         aot_value = lut( za_rad, az_rad )
         sensitivity[ i , j ] = aot_value
         

   import fits_beam
   fits_beam.save_fits( sensitivity , out_fitsname )
         

 

 
def parse_options(idx):
   usage="Usage: %prog [options]\n"
   usage += "Different plots / data can be obtainted:\n"
   usage += " 1/ to plot AoT vs. freq. for a given pointing direction (az,za) and lst time use options --azim_deg, --za_deg, --lst_hours\n"
   usage += " 2/ to plot AoT vs. time for a specified frequency use options --freq_mhz, --lst_start, --lst_end\n"
   usage += " 3/ to plot sensitivity map over the entire sky at a specified LST and frequency specify options : --lst_hours and --freq_mhz\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-p','--plot',action="store_true",dest="do_plot",default=False, help="Plot")
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


   # output file :
   parser.add_option('-o','--out_file','--outfile','--outout_file',dest="output_file",default=None, help="Full path to output text file basename (X or Y is added at the end) [default %default]" )
 

   (options, args) = parser.parse_args(sys.argv[idx:])
   
   print "###############################################################"
   print "PARAMATERS : "
   print "###############################################################"
   print "Do plotting                = %s" % (options.do_plot)
   print "Pointing direction (az,za) = (%.4f,%.4f) [deg]" % (options.azim_deg,options.za_deg)
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
      print "Unix time range            = %.2f - %.2f" % (options.unixtime_start,options.unixtime_end)
   else :
      print "Unix time range not specified"
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
    elif options.lst_hours is not None and options.frequency_mhz is not none :
       # plotting sensitivity map of the whole sky for a given frequnecy and pointing direction :
       (azim_x,za_x,aot_x,sefd_x,azim_y,za_y,aot_y,sefd_y) = sensitivity_db.get_sensitivity_map( options.frequency_mhz, options.lst_hours )
       
       if ( azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None ) or ( azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None ) :
#           if options.output_file is not None :
#              if azim_x is not None and za_x is not None and aot_x is not None and sefd_x is not None :
                  # save map of sensitivity to file 
              
#              if azim_y is not None and za_y is not None and aot_x is not None and sefd_y is not None :
                  # save map of sensitivity to file
                  
           # test :
           if options.do_plot :
              plot_sensitivity_map( azim_x, za_x, aot_x )

    # TODO :
    # option to plot map of the sky in zenith projection showing sensitivity across the whole sky
    #    help : https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
    #           https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # 
    # option to plot sensitivity at a given pointing direction and frequency and over a specified time range 
    
        
    