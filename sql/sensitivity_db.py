import sqlite3
from optparse import OptionParser

# plotting :
from pylab import *
import numpy
import matplotlib.pyplot as plt


# script for quering SQLITE3 or PostgreSQL databases for sensitivity values :
# HELP : python : https://www.sqlitetutorial.net/sqlite-python/ , https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/


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
                             station="EDA2", db_base_name="ska_station_sensitivity", db_path="./", 
                             db_lst_resolution=0.5, db_ang_res_deg=5.00 ) :
                             
    dbname_file = "%s/%s_%s.db" % (db_path,db_base_name,station)   
    
    conn = create_connection_sqlite3( dbname_file )


    cur = conn.cursor()
    szSQL = "SELECT id,azim_deg,za_deg,frequency_mhz,polarisation,lst,unixtime,gpstime,sensitivity,t_sys,a_eff,t_rcv,t_ant,array_type,timestamp,creator,code_version FROM Sensitivity WHERE ABS(lst-%.4f)<%.4f AND (za_deg-%.4f)<%.4f AND (az_deg-%.4f)<%.4f" %  (lst_hours,db_lst_resolution, za_deg, db_ang_res_deg, az_deg, db_ang_res_deg )
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

def plot_sensitivity( freq_x, aot_x, freq_y, aot_y, point_x='go', point_y='rx' ):

   plt.figure()
   plt.plot( freq_x, aot_x, point_x )
   plt.plot( freq_y, aot_y, point_y )
   plt.xlabel('Frequency (MHz)')
   plt.ylabel('Sensitivity A/T [m^2/K]')
   plt.grid()
   plt.show()

 
def parse_options(idx):
   usage="Usage: %prog [options]\n"
   parser = OptionParser(usage=usage,version=1.00)

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-p','--plot',action="store_true",dest="do_plot",default=False, help="Plot")
   parser.add_option('-a','--azim','--azim_deg',dest="azim_deg",default=None, help="Pointing direction azimuth in degrees [default %default]",metavar="float",type="float")
   parser.add_option('-z','--za','--za_deg',dest="za_deg",default=None, help="Pointing direction zenith distance in degrees [default %default]",metavar="float",type="float")
   parser.add_option('-l','--lst','--lst_hours',dest="lst_hours",default=None, help="Local sidreal time in hours [default %default]",metavar="float",type="float")

   # output file :
   parser.add_option('-f','--out_file','--outfile','--outout_file',dest="output_file",default=None, help="Full path to output text file basename (X or Y is added at the end) [default %default]" )
 

   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)

def save_output_file( freq, aot, pol, out_file_base ) :
   outfile_pol = ( "%s_%s.txt" % (out_file_base,pol) ) 
   out_f = open( outfile_pol , "w" )
   
   len = freq.shape[0]
   for i in range(0,len) :
      line = "%.4f %.8f\n" % (freq[i],aot[i])
      
      out_f.write( line )
   
   out_f.close()
 
if __name__ == "__main__":
    (options, args) = parse_options(1)

    if options.azim_deg is not None and options.za_deg is not None and options.lst_hours is not None :
       (out_freq_x,out_aot_x,out_sefd_x,out_freq_y,out_aot_y,out_sefd_y) = get_sensitivity_azzalst( options.azim_deg, options.za_deg, options.lst_hours )

       if options.output_file is not None :
          save_output_file( out_freq_x,out_aot_x, "X" , options.output_file )
          save_output_file( out_freq_y,out_aot_y, "Y" , options.output_file )

       if options.do_plot :
          plot_sensitivity( out_freq_x,out_aot_x, out_freq_y,out_aot_y )
          # do plot AoT vs. Freq :
        
    