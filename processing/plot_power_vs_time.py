#!/usr/bin/python

# plot power vs. time , based on ~/aavs-calibration/monitoring/plot_delay.py

# import pyfits
from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except ImportError:
    print("ERROR : could not import astropy.io.fits -> trying pyfits")
    import pyfits

# Test if X-windows / DISPLAY is not needed :
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import pylab
import math as m
from array import *
import matplotlib.pyplot as plt
import numpy as np
import string
import sys
import os
import errno
from optparse import OptionParser,OptionGroup
from pylab import *
from datetime import datetime
# from fit_line import fit_hor_line

# import rcv
# import metadata_auto

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      pass
            

def parse_options(idx):
   usage="Usage: %prog [options]\n"
   usage+='\tPlots visibilities for given tile and all other tiles\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-p','--plot','--show','--gui',action="store_true",dest="do_plot",default=False, help="Plot data and fit")
   parser.add_option('-i','--iter','--n_iter',dest="n_iter",default=1, help="Number of iteration to remove outliers from the fit",metavar="int",type="int")
   parser.add_option('-m','--max_resid','--limit',dest="limit",default=10, help="Fitting until residuals below this limit",metavar="float",type="float")
   parser.add_option('-o','--outfile','--out',dest="outfile",default=None, help="Name of output file [default %default]")
   parser.add_option('-s','--sqlfile','--sql',dest="sqlfile",default=None, help="Name of output SQL file [default %default]")
   parser.add_option('-c','--comment',dest="comment",default=None, help="Comment to be put in output file [default %default]")
   parser.add_option('--no_png',action="store_false",dest="save_png",default=True, help="Save png file")
   parser.add_option('-a','--amp_dir',dest="amp_dir",default="../amp/", help="Directory where gain txt files are stored [default %default]")
   parser.add_option('--obsid',dest="obsid",default=-1, help="Observations ID for SQL INSERT",metavar="int",type="int")
   parser.add_option('--tileid',dest="tileid",default=-1, help="Tile ID",metavar="int",type="int")
   parser.add_option('--last',dest="mean_last_n",default=4, help="Show mean of last N delays [default %default]",metavar="int",type="int")
   parser.add_option('--metafits',dest="metafits",default=None, help="Metafits file [default %default]")
   parser.add_option('--outdir','--dir',dest="outdir",default="images/", help="Output directory [default %default]")
      
   parser.add_option('-u','--uxtime','--unixtime','--unix_time',action="store_true",dest="unixtime",default=False, help="Interprete time column as unix time [default %s]")

   parser.add_option('--x_axis_title',dest="x_axis_title",default="Local time", help="X-axis title [default %default]")
   parser.add_option('--y_axis_title',dest="y_axis_title",default=None, help="Y-axis title [default %default]")
   parser.add_option('--power_unit',dest="power_unit",default="[?]", help="Delay units [default %default]")

   parser.add_option('--y_min',dest="y_min",default=None, help="Lower limit on the Y-axis [default %default]",type="float")
   parser.add_option('--y_max',dest="y_max",default=None, help="Upper limit on the Y-axis [default %default]",type="float")
   parser.add_option('--y_auto_median',dest="y_auto_median",default=False,action="store_true", help="Plot around median use --y_auto_median_range to specify the limits around median [default %]")
   parser.add_option('--y_auto_median_range',dest="y_auto_median_range",default=1.00,help="Plot around median range use --y_auto_median to enable [default %]")
   parser.add_option('--multiplier',dest="multiplier",default=None, help="Delay multiplier for example 1000.00 to convert from us -> nanosecond [default %default]",type="float")

   # --mean_stddev
   parser.add_option('--mean_stddev',dest="mean_stddev",default=False,action="store_true", help="Plot mean/stddev delay in txt file [default %]")
   
   # publication version :
   parser.add_option('--publication',dest="publication",default=False,action="store_true", help="Publication quality [default %s]")

   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)


def plotfile( filename,
              outdir      = "images/",
              do_gui      = True,
              mean_last_n = 4,
              metafits    = None,               
              unixtime = False,
              x_axis_title = 'Local time',
              y_axis_title = None,
              power_unit   = '[?]',
              y_min_ext    = None,
              y_max_ext    = None,
              multiplier   = None,
              comment      = None,
              y_auto_median = False,
              y_auto_median_range = 1.00,
              publication = False,
              db          = False
            ) :
              
   if y_axis_title is None :
       y_axis_title = ( "Power %s" % power_unit ) 
                     
   pngfile=filename.replace('.txt', '.png' )
   tile_id=filename.replace('.txt','')
   alldata = np.loadtxt(filename,usecols=[0,1])

   if len(alldata)>0 :
      uxtime=alldata[:,0]   
#      GPSfromUTC = (datetime(1980,1,6) - datetime(1970,1,1)).total_seconds()
#      uxtime=gpstime+GPSfromUTC
      uttime=[datetime(1980,1,6)]*uxtime.size
#      if unixtime :
#         uxtime=gpstime

      for i in range(0,len(uxtime)) :
         uttime[i]=datetime.utcfromtimestamp(uxtime[i])
    
      power=alldata[:,1]
      if db :
         power = 10.00*numpy.log10(power)
      
      y_min=min(power)-2
      y_max=max(power)+2

      x_min=min(uttime)
      x_max=max(uttime)

      fig=plt.figure()
      ax = fig.add_subplot(1,1,1) # create axes within the figure : we have 1 plot in X direction, 1 plot in Y direction and we select plot 1
      plt.ylim([y_min,y_max])
#      line_x, = plt.plot(uttime,x_delay_m, linestyle='None', marker='x', color='blue', markersize=10, label='X pol.')
#      line_y, = plt.plot(uttime,y_delay_m, linestyle='None', marker='x', color='red', markersize=10, label='Y pol.')
      if publication :
         line_x, = plt.plot(uttime,power, linestyle='None', marker='x', color='black', markersize=5, label='Power')
      else :
         line_x, = plt.plot(uttime,power, linestyle='None', marker='x', color='black', markersize=5, label='Power')
      plt.gcf().autofmt_xdate()
   
      ax.set_xlabel( x_axis_title )
      ax.set_ylabel( y_axis_title ) # r to treat it as raw string 
      
      title = filename
      if comment is not None :
         title = title + " , "
         title = title + comment 
         
      if publication :
         # Just antenna :
         title = comment   
         
      ax.set_title( title )
   


#      (B_x) = fit_hor_line(gpstime,x_delay_m,save_png=False,do_plot=False,limit_in_sigma=5,n_iter=10)
#      print("Fitted delay X pol = %.2f" % (B_x))

      # X-pol :
#      desc_x="Fitted X-pol delay = %.2f %s" % (B_x,power_unit)
#      if not publication :
#         plt.text((x_min), y_max-(y_max-y_min)*0.05, desc_x, fontsize=15, color='blue')   
#      x_delay_last=x_delay_m[-1]
#      if len(x_delay_m) < mean_last_n :
#         mean_last_n = len(x_delay_m)
      
   else :
      # no data case 
      fig=plt.figure()
      ax = fig.add_subplot(1,1,1) # create axes within the figure : we have 1 plot in X direction, 1 plot in Y direction and we select plot 1
      plt.ylim([0,1])
      x_min = datetime.utcfromtimestamp(1499050634-1*86400)
      x_max = datetime.utcfromtimestamp(1499050634)
      x_center = datetime.utcfromtimestamp(1499050634-86400)
      plt.xlim([x_min,x_max])
      ax.set_xlabel( x_axis_title )
      ax.set_ylabel( y_axis_title ) # r to treat it as raw string 
      ax.set_title(filename)
      plt.text(x_center, 0, "ERROR - no data from file: \n      " + filename, fontsize=30, color='red')
   
   

   mkdir_p(outdir)
   pngfile=outdir + "/" + pngfile
   plt.savefig(pngfile)   
   print("Saved file %s" % (pngfile))
   
   if do_gui :
      print("showed ?")
      plt.show()
      print("showed ?")
      
def plot_mean_stddev( filename,
              outdir      = "images/",
              do_gui      = True,
              mean_last_n = 4,
              metafits    = None,               
              unixtime = False,
              x_axis_title = 'Antenna Index',
              y_axis_title = None,
              power_unit   = '[?]',
              y_min_ext    = None,
              y_max_ext    = None,
              multiplier   = None,
              comment      = None,
              y_auto_median = False,
              y_auto_median_range = 1.00
            ) :
              
#   plt.rc('axes', titlesize=40)
#   plt.rc('axes', labelsize=40)
#   plt.grid(True)
   
                 
   if y_axis_title is None :
       y_axis_title = ( "Delay %s" % power_unit ) 
                     
   pngfile=filename.replace('.txt', '.png' )
   tile_id=filename.replace('.txt','')
   alldata = np.loadtxt(filename,usecols=[0,1,2,3,4])
   rcv_list = None
   if metafits is not None :
      print("Metafits = %s -> reading RCV-Tile mapping" % (metafits))
      rcv_list = metadata_auto.get_tile_rcv_mapping( metafits )

   if len(alldata)>0 :
      if len(alldata.shape) <= 1 :
         alldata = alldata.reshape( (1,3) )
         
      ants = alldata[:,0]   
      mean_x_delay_ns = alldata[:,1]
      stddev_x_delay_ns = alldata[:,2]
      mean_y_delay_ns = alldata[:,3]
      stddev_y_delay_ns = alldata[:,4]
      
      if multiplier is not None and multiplier > 0 :
          mean_x_delay_ns = mean_x_delay_ns * multiplier
          stddev_x_delay_ns = stddev_x_delay_ns * multiplier
          mean_y_delay_ns = mean_y_delay_ns * multiplier
          stddev_y_delay_ns = stddev_y_delay_ns * multiplier
   
      y_min=min(min(mean_x_delay_ns),min(mean_x_delay_ns))-2
      y_max=max(max(mean_y_delay_ns),max(mean_y_delay_ns))+2
      print("y_min = %.4f, y_max = %.4f" % (y_min,y_max))

      # overwrite to skip extreme points (at least for testing time to adjust phases to 0)
      if y_min < -50 :
         y_min = -50
      if y_max > 50 :
         y_max = 50
         
      if y_min_ext is not None and y_max_ext is not None :
          y_min = y_min_ext 
          y_max = y_max_ext

      if y_auto_median :
          y_median = np.median( y_delay_m )
          y_min = y_median - y_auto_median_range
          y_max = y_median + y_auto_median_range
          print("Using median = %.2f to plot in range %.3f - %.3f" % (y_median,y_min,y_max))

   
      x_min=1
      x_max=len( ants )

#      fig=plt.figure()
      fig = plt.figure(figsize=(60,20)) 
      ax = fig.add_subplot(1,1,1) # create axes within the figure : we have 1 plot in X direction, 1 plot in Y direction and we select plot 1
      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
         item.set_fontsize(30)
      
      plt.ylim([y_min,y_max])
      print("plt.ylim([%.4f,%.4f])" % (y_min,y_max))
#      line_x, = plt.plot(ants,mean_x_delay_ns, linestyle='None', marker='o', color='blue', markersize=3, label='X pol.')
#      line_y, = plt.plot(ants,mean_y_delay_ns, linestyle='None', marker='+', color='red', markersize=3, label='Y pol.')
      # line_x, = 
      plt.errorbar( ants, mean_x_delay_ns, yerr=stddev_x_delay_ns, xerr=None, linestyle='None', marker='o', color='blue' , markersize=5, label='X pol.')
      # line_y, = 
      plt.errorbar( ants, mean_y_delay_ns, yerr=stddev_y_delay_ns, xerr=None, linestyle='None', marker='+', color='red'  , markersize=5, label='Y pol.')

#      plt.legend(bbox_to_anchor=(0.85, 0.95),loc=3,handles=[line_x, line_y])
#      plt.legend([line_x, line_y], ['X pol.', 'Y pol.'])
#      plt.legend(handles=[line_x])
      # beautify the x-labels
      plt.gcf().autofmt_xdate()
      
#  matplotlib.pyplot.errorbar(x, y, yerr=None, xerr=None, fmt='', ecolor=None, elinewidth=None, capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, *, data=None, **kwargs)[source]
   
      ax.set_xlabel( x_axis_title , fontsize=30 )
      ax.set_ylabel( y_axis_title , fontsize=30 ) # r to treat it as raw string 
      
      title = filename
      if comment is not None :
         title = title + " , "
         title = title + comment          
         
      ax.set_title( title )
      plt.grid(True,which='both')
   


   else :
      # no data case 
      fig=plt.figure()
      ax = fig.add_subplot(1,1,1) # create axes within the figure : we have 1 plot in X direction, 1 plot in Y direction and we select plot 1
      plt.ylim([-180,180])
      x_min = datetime.utcfromtimestamp(1499050634-1*86400)
      x_max = datetime.utcfromtimestamp(1499050634)
      x_center = datetime.utcfromtimestamp(1499050634-86400)
      plt.xlim([x_min,x_max])
      ax.set_xlabel( x_axis_title )
      ax.set_ylabel( y_axis_title ) # r to treat it as raw string 
      ax.set_title(filename)
      plt.text(x_center, 0, "ERROR - no data from tile: \n      " + filename, fontsize=30, color='red')
   
   

   mkdir_p(outdir)
   pngfile=outdir + "/" + pngfile
   plt.savefig( pngfile , dpi=80 )   
   print("Saved file %s" % (pngfile))

   
   if do_gui :
      plt.show()

   
#   plt.plotfile(filename,(0,1),delimiter=" ",names="Frequency [MHz],T[K]",newfig=False)
#   plt.savefig(pngfile)
   

def main() :
   # PARSE COMMAND LINE :
   filename="EDA2_power_vs_time.txt"
   if len(sys.argv) > 1:
      filename = sys.argv[1]
   (options, args) = parse_options(1)      

   plotfile( filename,
                do_gui=options.do_plot,
                mean_last_n=options.mean_last_n,
                metafits=options.metafits,
                unixtime     = options.unixtime,
                x_axis_title = options.x_axis_title,
                y_axis_title = options.y_axis_title,
                power_unit   = options.power_unit,
                y_min_ext    = options.y_min,
                y_max_ext    = options.y_max,
                multiplier   = options.multiplier,
                comment      = options.comment,
                outdir       = options.outdir,
                y_auto_median = options.y_auto_median,
                y_auto_median_range = options.y_auto_median_range,
                publication = options.publication
              )
#   else :
#      print "ERROR : file %s is empty !" % filename

if __name__ == "__main__":
   main()
