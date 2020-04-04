import h5py
import sys
import time
from optparse import OptionParser,OptionGroup

def parse_options(idx):
   usage="Usage: %prog hdf5fits_station_beam.py station_beam.hdf5 [options]\n"
   usage+='\tDump station beam power to a text file\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-c','--channel','--ch',dest="channel",default=4, help="Channel to dump (value <0 - means mean of all channels) [default %default]",metavar="int",type="int")
   parser.add_option('-p','--pol',dest="polarisation",default=0, help="Polarisation [default %default]",metavar="int",type="int")
   parser.add_option('-o','--outfilebase','--out_file_base','--out_file_basename',dest="out_file_basename",default="power_vs_time_ch%d_%s.txt", help="Output file name template [default %default]" )
   parser.add_option('-l','--last_n_seconds','--interval',dest="last_n_seconds",default=630720000,help="Save onle the last N seconds [default %default ~= infinitly in the past, i.e. all data]",type="int" )
   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)



hdf5file="stationbeam_integ_0_20190326_84842_0.hdf5"
if len(sys.argv) > 1: 
  hdf5file = sys.argv[1]

(options, args) = parse_options(1)

channel=options.channel
pol_name="X"
if options.polarisation == 1 :
   pol_name="Y"

f = h5py.File( hdf5file )

print "keys = %s" % (f.name)
print "keys = %s" % f.keys()


pol_string=( "/polarization_%d" % (options.polarisation) )
data=f[pol_string]['data']
times=f['/sample_timestamps']['data']

print "Data.shape = %d x %d" % (data.shape[0],data.shape[1])

n_timesteps=times.shape[0]

current_uxtime = time.time()
print "Plotting power after uxtime = %d" % (current_uxtime - options.last_n_seconds)

outfile_name = options.out_file_basename % (channel,pol_name)
out_f = open( outfile_name , "w" )
for t in range(0,n_timesteps) :
   # mean = ( data[t][0] + data[t][1] + data[t][2]  + data[t][3]  + data[t][4]  + data[t][5]  + data[t][6]  + data[t][7] ) / 8
   
   # use 4th channel 204 :
   mean = 0.00
   if channel >=0 :
       mean = data[t][channel] 
   else :
       # if channel not specified ( < 0 ) -> use mean of 8 channels :
       mean = ( data[t][0] + data[t][1] + data[t][2]  + data[t][3]  + data[t][4]  + data[t][5]  + data[t][6]  + data[t][7] ) / 8
       
   if mean > 0 :    
       if times[t][0] > ( current_uxtime - options.last_n_seconds) :
           line = "%.4f %.4f\n" % (times[t][0],mean)
           out_f.write( line )
   else :
       print "mean = %.4f <= 0 -> skipped" % (mean)
   

out_f.close()
f.close()
