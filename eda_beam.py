#!/usr/bin/python
"""Tools for calculating the complex voltage response and power pattern of the EDA with MWA dipoles
Randall Wayth & Kim Steele. June 2016.
"""
from __future__ import print_function
# import pdb

import numpy,logging
import astropy.io.fits as pyfits
import os,sys
from scipy import interpolate
import math
import fits_beam

vel_light=2.99792e8
DQ=435e-12*vel_light  # delay quantum is distance light travels for 1 quantum (meters). This for MWA beamformers
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger(__name__)   # default logger level is WARNING
filesdir=os.path.dirname(__file__)
if (len(filesdir)==0):
    filesdir='.'
logger.debug("Looking for data files in "+filesdir)

# OLD locations as RW/KS 
# xs=[-14.78,10.16,-2.13,12.07,4.96,-2.43,5.36,0.91,4.07,-14.29,2.40,0.30,0.87,0.51,7.81,-12.10,-1.08,5.16,-14.15,12.60,0.59,-0.51,1.42,-2.32,-13.37,-9.27,-2.05,-4.14,-5.84,2.90,13.77,4.39,3.64,-5.31,-9.17,-1.43,4.35,8.57,11.07,-13.05,-0.96,-6.89,5.68,-6.46,15.72,15.03,-2.79,14.22,-6.81,-16.33,-5.39,-15.85,6.61,2.92,14.61,6.81,-16.10,-2.64,-4.83,16.43,-8.33,9.15,-15.52,-14.88,-12.58,-0.84,6.77,0.08,6.43,0.44,-9.23,-8.18,-3.67,9.99,-6.01,5.29,-2.91,-5.50,-4.75,1.26,2.76,-13.63,13.58,-8.78,-12.51,2.18,3.08,6.12,13.74,-7.80,9.71,11.36,17.52,12.56,-10.05,16.81,3.63,14.17,7.44,-10.85,-11.24,1.08,-4.73,-12.33,-0.32,-8.35,2.50,7.84,16.63,-8.89,2.84,12.55,-8.59,-10.79,-5.77,-6.40,-7.44,-3.47,-4.48,-0.76,12.17,4.10,2.95,-8.21,1.90,-6.57,6.00,-2.83,-11.01,9.70,-4.39,14.71,-10.54,-6.81,11.52,8.17,-1.03,-3.99,13.88,2.83,12.32,-10.64,-13.05,-6.42,9.92,0.12,8.41,10.97,-6.72,-9.25,1.64,-0.14,8.20,5.08,2.42,-4.48,11.50,-0.30,-16.94,8.82,11.25,12.67,-2.81,8.45,5.73,-11.71,8.66,3.30,2.53,-16.08,-12.25,6.79,16.21,-4.14,4.64,0.06,-11.90,0.99,9.33,-0.36,4.52,-9.91,-15.71,2.61,4.23,-7.88,-2.22,5.65,-1.45,13.82,9.90,-2.15,-14.20,-12.95,11.30,11.11,-2.87,2.67,4.93,-5.49,-11.92,4.34,7.32,12.81,-14.92,1.73,-10.24,-3.93,-17.33,4.78,6.94,-12.38,-3.19,-9.75,10.25,1.42,1.31,-2.13,-8.37,-1.54,-15.90,11.17,-11.19,-7.03,-10.17,-15.26,-14.08,-14.20,-3.21,-16.99,-12.89,6.97,-9.43,7.73,6.18,0.00,-4.78,1.87,10.65,-10.85,4.95,-7.48,14.08,-9.31,-3.08,14.70,8.28,-7.56,-2.49,-9.02,9.11,10.46,9.44,-15.16,-13.92,10.11]
# ys=[-7.66,1.94,5.53,11.94,-1.83,-3.17,1.08,6.09,-13.26,6.46,-5.87,-11.44,-2.59,-15.74,4.89,8.11,6.47,-9.88,-5.87,0.12,9.51,8.20,11.83,-14.02,11.33,0.91,10.67,-13.98,0.33,-14.42,3.27,-0.13,4.66,8.11,-1.82,-10.00,3.32,13.83,-3.63,4.62,-6.51,-5.17,-15.62,-15.20,4.80,-1.85,-15.64,-7.27,-11.83,-5.91,10.71,5.17,-2.73,12.93,8.21,-12.77,-3.13,8.31,-11.45,-3.66,6.29,-0.45,1.31,9.44,-8.74,9.55,7.42,-4.27,4.95,1.70,10.28,8.12,14.28,10.66,-9.02,-5.37,-5.78,-12.84,16.03,-17.55,-8.86,0.62,-9.46,-10.82,9.76,-0.66,10.18,-8.24,4.85,-14.47,-9.39,8.98,-1.12,-4.09,-3.64,0.94,-11.94,-5.11,-14.78,6.04,11.32,2.96,-4.81,-12.51,16.11,11.72,15.34,-9.88,3.04,-7.65,6.38,-7.21,3.03,-9.32,-3.31,6.78,-0.80,-12.70,9.38,13.75,2.27,14.31,-4.12,-4.08,-11.24,9.32,12.15,3.64,1.06,-14.63,2.65,-3.58,2.94,13.83,0.99,1.06,-1.50,5.52,0.83,-17.29,-11.07,-6.96,-4.39,-1.86,-7.79,-7.78,-8.35,-6.51,2.53,4.47,16.94,10.92,-4.61,8.31,-2.43,-9.09,13.20,-0.25,-0.31,9.27,6.69,-1.84,1.01,-12.51,14.67,-10.74,12.27,8.60,1.64,6.98,12.35,1.28,-6.32,-16.38,-16.87,-13.51,-5.50,-5.50,-2.73,3.26,-8.23,8.89,2.74,-7.41,5.97,15.02,15.31,6.40,-17.04,10.85,7.15,13.18,-3.40,2.78,4.58,-9.68,17.40,-15.92,-6.80,12.21,3.82,17.07,-5.72,8.27,3.98,4.27,-11.68,-7.68,3.12,10.96,-1.04,-0.74,-1.36,12.97,3.53,-13.92,13.61,-7.71,-5.61,-4.80,-1.77,-2.13,-2.17,5.40,-13.67,-0.42,-1.42,-8.95,12.10,1.31,-7.21,16.05,7.31,-11.31,9.99,12.54,-0.81,8.99,-0.77,9.95,-11.43,-9.75,-0.64,-9.28,-11.05,6.20,6.59,10.39,-9.09,-12.89,5.34,12.21,-10.84,-4.62,-10.65,-12.47]

# NEW LOCATIONS FROM locations.txt file in MandC (bigdas):
xs=[-2.862,-2.183,-4.774,-0.325,-3.654,-0.74,-2.133,1.642,-3.187,-6.792,-7.917,0.004,1.293,2.504,-5.492,-2.05,-14.874,-12.507,-13.352,-16.091,-14.296,-12.116,-12.253,-10.829,-11.253,-15.856,-9.905,-13.061,-10.942,-14.912,-9.254,-9.453,4.366,4.112,5.753,6.93,2.913,5.993,8.561,1.41,4.804,8.656,3.087,6.197,-0.156,10.45,11.513,0.591,-17.351,-15.698,-17.001,-15.497,-16.951,-15.246,-12.946,-13.669,-15.917,-11.908,-14.1,-12.389,-16.111,-11.003,-10.547,-14.209,1.263,2.857,0.514,2.636,-1.431,4.655,2.904,1.392,0.043,-2.793,5.717,-2.289,4.037,-4.145,3.636,0.298,13.789,12.036,14.608,12.816,11.344,10.01,14.716,11.264,8.784,9.881,13.764,15.705,11.309,8.304,9.085,13.758,-12.376,-11.71,-10.235,-13.916,-10.167,-8.988,-10.825,-12.589,-8.816,-14.217,-9.304,-7.816,-12.915,-14.807,-6.786,-7.476,17.52,16.803,15.063,16.444,14.114,14.692,13.904,16.656,12.648,12.608,14.149,16.192,12.552,12.154,11.491,11.141,9.697,10.084,7.468,8.457,6.789,9.453,7.736,12.308,7.844,11.088,9.7,4.963,8.404,13.594,5.147,9.904,-16.362,-15.17,-14.164,-13.036,-11.915,-10.659,-11.216,-10.031,-8.914,-8.337,-9.162,-8.202,-6.899,-9.261,-7.447,-6.411,-6.433,-5.486,-4.14,-3.463,-4.842,-3.089,-5.988,-4.491,-1.428,-2.497,-3.959,-2.152,1.901,0.101,-2.897,-0.963,-9.756,-8.358,-7.532,-6.541,-5.384,-8.199,-4.489,-5.288,-8.358,-6.388,-7.051,-2.643,-9.245,-3.992,-0.854,-8.585,14.218,12.549,10.947,11.045,8.206,9.341,7.297,10.646,6.087,9.145,6.595,5.279,4.938,6.971,4.492,10.143,10.252,7.792,8.143,6.43,6.812,6.752,5.672,5.335,4.31,4.26,3.671,4.404,5.095,4.962,2.848,2.54,2.744,2.604,2.38,0.99,2.957,0.1,-1.545,2.396,0.84,-2.469,2.194,-1.058,-4.727,-0.3,-3.199,-5.754,-2.148,0.905,-1.105,-0.502,-5.833,0.467,1.105,-2.835,-4.412,-6.684,-2.8,3.311,-0.355,1.755,-4.789,1.85]
ys=[17.403,15.282,16.063,16.093,14.285,13.725,13.182,16.927,12.103,13.849,14.994,12.543,13.587,15.328,12.189,10.651,9.453,9.77,11.355,6.981,6.441,8.103,12.363,9.952,11.288,5.157,8.9,4.605,6.063,3.987,10.288,7.297,17.057,14.292,14.638,16.054,12.931,12.157,13.846,11.835,10.952,12.258,10.217,9.993,10.907,12.217,13.2,9.513,3.087,2.733,1.286,1.285,-0.293,-0.414,2.799,0.581,-1.758,3.801,-1.397,-0.759,-3.168,1.05,2.946,-3.415,-17.552,-17.29,-15.749,-15.901,-17.048,-16.851,-14.392,-13.896,-13.484,-15.662,-15.612,-14.021,-13.253,-16.397,-11.936,-11.452,10.852,11.949,8.211,8.259,9.023,10.661,6.196,6.701,9.238,7.153,4.841,4.789,4.596,6.612,5.353,3.244,-12.504,-10.743,-11.699,-10.641,-13.662,-12.891,-9.302,-8.739,-10.801,-8.943,-9.294,-14.459,-7.2,-7.635,-11.855,-9.758,-1.107,0.954,-1.841,-3.658,-0.647,-3.579,0.828,3.05,-1.83,0.091,-5.088,-6.311,-4.083,2.245,1.015,-2.168,-14.654,-12.449,-14.803,-12.485,-12.76,-10.85,-11.288,-11.036,-9.897,-9.689,-9.414,-11.44,-8.345,-9.471,-9.886,-7.818,-5.906,-4.606,-5.852,-4.385,-5.518,-6.96,-2.176,-3.662,-7.639,-5.609,-1.814,-4.096,-5.163,0.886,-0.821,-1.868,-15.194,-12.852,-13.982,-12.712,-11.461,-11.059,-8.997,-9.096,-10.007,-9.105,-7.707,-7.71,-11.239,-7.803,-5.801,-6.501,12.952,11.699,10.397,9.288,10.707,8.095,9.392,8.101,6.294,6.794,5.4,8.304,4.442,5.494,9.531,3.049,-7.234,-7.195,-6.484,-3.639,-4.626,-2.754,-5.688,-0.766,-8.256,-0.447,-2.761,-5.354,-6.785,-1.048,-8.252,1.95,3.551,4.903,1.041,4.962,1.308,7.399,6.397,1.085,3.295,5.955,4.645,-0.144,8.319,-1.851,6.386,1.642,-8.864,-7.401,-5.842,-5.501,-4.083,-4.243,-4.8,-2.471,-2.595,-3.145,-0.649,-1.511,-4.804,-0.25,-1.347,-3.326,5.556,6.103,6.469,8.217,0.338,1.702,2.962,3.644,2.641,2.548,1.001,8.612,3.254,4.251,-0.795,9]
n_dipoles_eda=256

def read_antenna_list( filename, 
                       start_column=1, 
                       overwrite=False # overwrite default EDA1 antenna list 
                     ):
   file = open( filename , 'r' )
   data = file.readlines()
   
   x_arr = []
   y_arr = []
   z_arr = []
   count = 0
   for line in data : 
       words = line.split() # was ' ' , but when none is provided -> all white-space characters are ok !
       print("DEBUG : words = %s (len = %d)" % (words,len(words)))

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

   print("Read %d / %d / %d of x, y, z positions from file %s" % (len(x_arr),len(y_arr),len(z_arr),filename))
   
   if overwrite :
      print("Overwritting %d antenna positions with a list from file %s" % (count,filename))
      global xs
      global ys
      global n_dipoles_eda 
      
      xs = numpy.copy( x_arr )
      ys = numpy.copy( y_arr )
      n_dipoles_eda = count     

   return (x_arr,y_arr,z_arr,count)



def distance( x, y, x0, y0 ) :
   return math.sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) )

# BF=4 (0-indexed) or 5 (1-indexed) excluded, based on /home/msok/Desktop/EDA/doc/EDA GP Mech analysis_1rev3.pdf
def is_broken_dipole( dip_x=-1.00, dip_y=-1.00 ) :
   return False

   if distance( dip_x, dip_y, 1.263345, -17.554271 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 2.832297,-17.291275 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 0.507340, -15.735098 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 2.670654, -15.919592 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, -1.450232, -17.040957 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 4.638593, -16.867783 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 2.897577, -14.417846 ) < 0.5 :
      return True
   
   if distance( dip_x, dip_y, 1.421085, -13.917899 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, 0.055905, -13.507449 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, -2.792250 , -15.638102 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, 5.684624, -15.620430 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, -2.315381, -14.016200 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, 4.071723, -13.255940 ) < 0.5 :
      return True
      
   if distance( dip_x, dip_y, -4.138596, -16.383160 ) < 0.5 :
      return True          
      
   if distance( dip_x, dip_y, 3.626534,-11.937362 ) < 0.5 :
      return True          
      
   if distance( dip_x, dip_y, 0.303865, -11.440366 ) < 0.5 :
      return True                 

      
   return False

class Dipole:
    """Provides a generic dual pol dipole object to support MWA beam models"""

    def __init__(self, type='short', height=0.29, len=0.74, lookup_filename=os.path.join(filesdir,'Jmatrix.fits'), gain=None,station_name="EDA", projection="zea" ):
        """General dipole object. Dual pol crossed dipole.
        Assumes a groundscreen with dipole height meters above ground.
        Supported types are 'short': analytic short dipole
                            'long':  analytic real length dipole
                            'lookup': lookup table
        For a short dipole, the length of the dipole is ignored.
        Gain is a 2x2 matrix with gain and direction independent crosstalk
        for the X (upper) and Y (lower) dipole resectively
        """
        assert type=='short' or type=='long' or type=='lookup' or type=='fits_beam', 'Unknown type %r' % type
        self.type=type
        self.height=height
        self.len=len
        self.station_name = station_name
        self.projection = projection
        if gain is None:
            self.gain = numpy.eye(2,dtype=numpy.complex64)
        else:
            self.gain=gain
        self.interp_freq = 0.0
        if type=='lookup':
            lookup=self.loadLookup(lookup_filename)

    def loadLookup(self,lookup_filename):
        """Load a dipole Jones response lookup table (FITS file)"""
        # data is a direct conversion of the output of the simulation so has redundant info.
        # we do all the conversion and stuff here
        # for reference, columns are:
        # theta phi  real(Jxt(t,p)) imag(Jxt(t,p)) real(Jxp(t,p)) imag(Jxp(t,p)) real(Jyt(t,p)) imag(Jyt(t,p)) real(Jyp(t,p)) imag(Jyp(t,p)))
        try:
            hdulist = pyfits.open(lookup_filename)
        except IOError:
            raise Exception('Cannot load Jones matrix file %s' % lookup_filename)
        nfreqs = len(hdulist)
        self.lookup_za = numpy.unique(hdulist[0].data[:,0]) # zenith angle
        self.lookup_ph = numpy.unique(hdulist[0].data[:,1]) # phi angle == 90 - az
        p = numpy.where
        nza = len(self.lookup_za)
        nph = len(self.lookup_ph)
        self.lookup = numpy.empty((nfreqs,nza,nph,2,2),dtype=numpy.complex64)
        freqs = []
        for i in range(nfreqs):
            hdu=hdulist[i]
            logger.debug('Loading J lookup matrix for freq '+str(hdu.header['FREQ']))
            freqs.append(hdu.header['FREQ'])
            jxt = hdu.data[:,2] + 1.0j*hdu.data[:,3]
            jxp = hdu.data[:,4] + 1.0j*hdu.data[:,5]
            jyt = hdu.data[:,6] + 1.0j*hdu.data[:,7]
            jyp = hdu.data[:,8] + 1.0j*hdu.data[:,9]
            self.lookup[i,:,:,0,0] = jxt.reshape((nph,nza)).transpose()
            self.lookup[i,:,:,0,1] = jxp.reshape((nph,nza)).transpose()
            self.lookup[i,:,:,1,0] = jyt.reshape((nph,nza)).transpose()
            self.lookup[i,:,:,1,1] = jyp.reshape((nph,nza)).transpose()
        logger.info('Loaded dipole Jones matrix lookup model from '+lookup_filename+' with '+str(nfreqs)+' freqs')
        self.freqs = numpy.array(freqs)
        logger.info('Supported frequencies (MHz): '+str(self.freqs/1e6))
        logger.debug("There are "+str(nza)+ " tabulated zenith angles: "+str(self.lookup_za))
        logger.debug("There are "+str(nph)+ " tabulated phi angles: "+str(self.lookup_ph))

    def getJones(self,az,za,freq,zenith_norm=True):
        """Return the Jones matrix for arrays of az/za for a given freq
        az and za are numpy arrays of equal length.
        Results are in corrds of az/za unit vectors such that
        output[j,i,:,:] is a 2x2 Jones matrix that maps
        za onto E-W     az onto E-W
        za onto N-S     az onto N-S
        By default, normalise to the zenith
        """
        if self.type=='short':
            return self.getJonesShortDipole(az,za,freq,zenith_norm)
        elif self.type=='lookup':
            return self.getJonesLookup(az,za,freq,zenith_norm)
        elif self.type=='fits_beam' :
            return self.getFitsBeam(az,za,freq)
        else:
            raise Exception("Dipole type %r is not implemented yet." % self.type)

    def getFitsBeam(self,az,za,freq) :
#        def get_fits_beam( azim_deg, za_deg, frequency_mhz, polarisation="X", station_name="EDA", simulation_path="$HOME/aavs-calibration/BeamModels/") :
        (beam_values_x,current_fits_beam_x,x_pixel_x,y_pixel_x) = fits_beam.get_fits_beam_multi( az, za, freq/1e6, "X", station_name=self.station_name, projection=self.projection )
        (beam_values_y,current_fits_beam_y,y_pixel_y,y_pixel_y) = fits_beam.get_fits_beam_multi( az, za, freq/1e6, "Y", station_name=self.station_name, projection=self.projection )
        
        print("DEBUG : fits_beam.get_fits_beam_multi returned beam_x = %.4f , beam_y = %.4f for pixel (%.1f,%.1f) at (az,za) = (%.4f,%.4f) - please verify" % (beam_values_x[0,0],beam_values_y[0,0],x_pixel_x[0,0],y_pixel_x[0,0],az[0,0],za[0,0]))

        result = numpy.empty((za.shape+(2,2)),dtype=numpy.complex64)
        
        # isotropic antenna :
        result[...,0,0] = beam_values_x.reshape( za.shape ) # this is jones so power beam need to be SQRT
        result[...,0,1] = 0.00 * beam_values_x.reshape( za.shape )
        result[...,1,0] = 0.00 * beam_values_x.reshape( za.shape )
        result[...,1,1] = beam_values_y.reshape( za.shape ) # this is jones so power beam need to be SQRT
        return result


    def getJonesLookup(self,az,za,freq,zenith_norm=True):
        """Return the Jones matrix for arrays of az/za for a given freq (Hz)
        this method interpolates from the tablulated numerical results loaded
        by the constructor"""
        # need to interpolate each of the 4 Jones elements separately and each
        # the real and imag separately (since interpolate.RectBivariateSpline)
        # apparently doesn't handle complex

        # find the nearest freq lookup table
        pos = numpy.argmin(numpy.abs(self.freqs - freq))
        if numpy.abs(self.freqs[pos] - freq) > 2e6:
            logger.warning("Nearest tabulated impedance matrix freq is more than 2 MHz away from desired freq.")
        logger.info("Selecting matrix for nearest freq "+str(self.freqs[pos]))

        # cache the interpolation functions
        if self.interp_freq != freq:
            self.interp_freq = freq
            logger.debug("Setting new cache lookup freq to "+str(self.freqs[pos]))
            self.i00_real = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,0,0].real)
            self.i00_imag = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,0,0].imag)
            self.i01_real = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,0,1].real)
            self.i01_imag = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,0,1].imag)
            self.i10_real = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,1,0].real)
            self.i10_imag = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,1,0].imag)
            self.i11_real = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,1,1].real)
            self.i11_imag = interpolate.RectBivariateSpline(self.lookup_za,self.lookup_ph,self.lookup[pos,:,:,1,1].imag)
            # determine normalisation factors. The simulations include all ph angles at za=0.
            # these are not redundant, and the ph value determines the unit vector directions of
            # both axes. We should normalise by where the result will be maximal.
            # For the E-W dipoles, the projection of the ZA unit vec will be max when 
            # pointing east, i.e. when ph=0. For the PH unit vec, this will be when ph=-90 or 90
            # For the N-S dipoles, projection of ZA onto N-S is max az ph=90 and
            # proj of ph onto N-S is max when ph=0

            # determine the indices of 0 and 90 degrees of phi in the tabulated values
            ph0 = numpy.where(self.lookup_ph == 0.0)
            ph90= numpy.where(self.lookup_ph == 90.0)
            th0 = numpy.where(self.lookup_za == 0.0)

            self.j00norm = self.lookup[pos,th0,ph0 ,0,0]
            self.j01norm = -self.lookup[pos,th0,ph90,0,1]   # use -90, not 90.
            self.j10norm = self.lookup[pos,th0,ph90,1,0]
            self.j11norm = self.lookup[pos,th0,ph0 ,1,1]

        ph_deg = 90.0 - az*180.0/numpy.pi   # ph in degrees
        za_deg = za*180.0/numpy.pi
        p = ph_deg < 0
        ph_deg[p] += 360.0
        j00 = self.i00_real.ev(za_deg.flatten(),ph_deg.flatten()) + 1.0j*self.i00_imag.ev(za_deg.flatten(),ph_deg.flatten())
        j01 = self.i01_real.ev(za_deg.flatten(),ph_deg.flatten()) + 1.0j*self.i01_imag.ev(za_deg.flatten(),ph_deg.flatten())
        j10 = self.i10_real.ev(za_deg.flatten(),ph_deg.flatten()) + 1.0j*self.i10_imag.ev(za_deg.flatten(),ph_deg.flatten())
        j11 = self.i11_real.ev(za_deg.flatten(),ph_deg.flatten()) + 1.0j*self.i11_imag.ev(za_deg.flatten(),ph_deg.flatten())

        result = numpy.empty((za.shape+(2,2)),dtype=numpy.complex64)
        result[...,0,0] = j00.reshape(za.shape)/self.j00norm
        result[...,0,1] = -j01.reshape(za.shape)/self.j01norm    # sign flip between az and phi
        result[...,1,0] = j10.reshape(za.shape)/self.j10norm
        result[...,1,1] = -j11.reshape(za.shape)/self.j11norm     # sign flip between az and phi
        return result

    def getJonesShortDipole(self,az,za,freq,zenith_norm=True):
        """Calculate the Jones matrix for a short dipole.
        This is defined by purely geometric projection of unit vectors
        on the sky onto the unit vector defined by the dipole's direction.
        """
        # check for input scalars vs arrays. If the input is scalar, then
        # convert it so that numpy array-based operations are possible
        if numpy.isscalar(az):
            logger.info("Converting scalar az input to a matrix")
            az = numpy.asarray(az,dtype=numpy.float32)
        if numpy.isscalar(za):
            logger.info("Converting scalar za input to a matrix")
            za = numpy.asarray(za,dtype=numpy.float32)
        assert az.shape == za.shape,"Input za and az arrays must have same dimension"

        # output array has 2x2 Jones matrix for every az/za point in input
        result = numpy.empty((za.shape+(2,2)),dtype=numpy.complex64)
        # apply the groundscreen factor, which is independent of az
        znorm=1.0
        gs = self.groundScreen(za,freq)
        if zenith_norm:
            logger.debug("Normalising response to zenith")
            gs /= self.groundScreen(0.0,freq)
        result[...,0,0] = numpy.cos(za)*numpy.sin(az)*gs
        result[...,0,1] = numpy.cos(az)*gs
        result[...,1,0] = numpy.cos(za)*numpy.cos(az)*gs
        result[...,1,1] = -numpy.sin(az)*gs
        return result

    def groundScreen(self,za,freq):
        """
        Calculate the groundscreen effect for an ideal infinite groundscreen
        given the dipole's height above the screen and the frequency (Hz)
        """
        l = vel_light/freq
        return numpy.sin(numpy.pi * (2.0*self.height/l) * numpy.cos(za))*2.0

    def __str__(self):
        return "Dipole. Type: "+self.type+". height: "+str(self.height)+"m. Gain: "+str(self.gain)


class ApertureArray:
    """Aperture array antenna object"""

    def __init__(self, dipoles=None, xpos=None, ypos=None, zpos=None):
        """Constructor for aperture array station. xpos and ypos are arrays with
        the coords of the dipoles (meters) in local coords relative to centre of
        the antenna looking down on the station. Ordering goes left to right, top to bottom,
        hence are offsets in east and north from the array phase centre.
        """
        global n_dipoles_eda
        
        if xpos is not None and ypos is not None :
            n_dipoles_eda = len( xpos )
            print("ApertureArray.__init__ set n_dipoles_eda = %d" % (n_dipoles_eda))
        
        if xpos is None :
           xpos = xs

        if ypos is None :
           ypos = ys                   
                   
        print("ApertureArray.__init__ initialising with positions :")
        print("\txs = %s" % (xpos))
        print("\tys = %s" % (ypos))      
        
        if dipoles == None:
            d = Dipole(type='short')
            dipoles = [d]*n_dipoles_eda
        self.dipoles = dipoles
        self.xpos = numpy.array(xpos)
        self.ypos = numpy.array(ypos)
        if zpos == None:
            self.zpos = self.xpos*0.0
        else:
            self.zpos = numpy.array(zpos)
            print("DEBUG : Initialised array with zpos != None : \n%s" % (self.zpos))

    def getPortCurrents(self,freq=155e6,delays=numpy.zeros((2,n_dipoles_eda),dtype=numpy.float32)):
        """
        Return the port currents on a tile given the freq (Hz) and delays (in units of the DQ)
        """
        n_dipoles=numpy.size(delays)/2
        lam = vel_light/freq
        phases = -2.0*numpy.pi*delays*(DQ/lam)
        ph_rot = numpy.cos(phases) + 1.0j*numpy.sin(phases)
        # this code ignores any dipole gain (and crosstalk) terms.
        # should FIXME it.
        #z_total=self.im.getImpedanceMatrix(freq) + numpy.eye(32)*self.lna_z.getZ(freq)
        #inv_z = numpy.linalg.inv(z_total)
        inv_z = numpy.eye(numpy.size(delays))/50.0  # hard-code 50 Ohm impedance for now
        port_current = numpy.dot(inv_z,ph_rot.reshape(numpy.size(delays))).reshape(2,n_dipoles)
        # apply electronic gains, including crosstalk
        # select out views of the big array and apply gains...
        # in port currents, Y dipoles are first, but in gains X dipoles are first. Hmmm.
        for i in range(n_dipoles):
            dipole=self.dipoles[i]
            temp = port_current[:,i]    # view into big array as [Y,X]. Changing temp changes the original...
            temp_r = numpy.dot(dipole.gain,temp[::-1])  # apply gains to get result as [X,Y]
            temp[::] = temp_r[::-1] # apply back to original as
        return port_current

    def getArrayFactor(self,az,za,freq=155e6,delays=numpy.zeros((2,n_dipoles_eda),dtype=numpy.float32)):
        """
        Get the scalar array factor response of the array for a given
        freq (Hz) and delay settings.
        az and za (radian) are numpy arrays of equal length defining a set
        of points to calculate the response for.
        delays is a 2D array of integer delay steps for the Y and X pol
        respectively.
        Result are in same coords as the az/za input arrays
        """
        lam = vel_light/freq
        port_current = self.getPortCurrents(freq,delays)

        # check for input scalars vs arrays. If the input is scalar, then
        # convert it so that numpy array-based operations are possible
        if numpy.isscalar(az):
            logger.info("Converting scalar az input to a matrix")
            az = numpy.matrix(az,dtype=numpy.float32)
        if numpy.isscalar(za):
            logger.info("Converting scalar za input to a matrix")
            za = numpy.matrix(za,dtype=numpy.float32)
        assert az.shape==za.shape, "Input az and za arrays must have same dimenions"

        # now calculate the array factor using these port currents. Note az and za
        # are likely to be arrays
        sz = numpy.sin(za)
        kx = (2.0*numpy.pi/lam)*numpy.sin(az)*sz
        ky = (2.0*numpy.pi/lam)*numpy.cos(az)*sz
        kz = (2.0*numpy.pi/lam)*numpy.cos(za)

        ax = numpy.zeros_like(az,dtype=numpy.complex64)
        ay = numpy.zeros_like(az,dtype=numpy.complex64)
        az = numpy.zeros_like(az,dtype=numpy.complex64)

        count_broken=0
        for i in range(len(self.xpos)):
            ph = kx*self.xpos[i] + ky*self.ypos[i] + kz*self.zpos[i]
                        
            # print "Antenna[%d] : phase = %.4f [deg]" % (i,ph*(180.00/numpy.pi))
            
            dipole_x = self.xpos[i]
            dipole_y = self.ypos[i] 
            # print "Dipole (x,y) = (%.2f,%.2f)" % (dipole_x,dipole_y)
            is_broken = is_broken_dipole( dip_x=dipole_x, dip_y=dipole_y )
            
            if is_broken :            
               count_broken = count_broken + 1 
            else :
               ax += port_current[1,i]*(numpy.cos(ph) + 1.0j*numpy.sin(ph))    # X dipoles
               ay += port_current[0,i]*(numpy.cos(ph) + 1.0j*numpy.sin(ph))    # Y dipoles
               
#               if ax.shape[0] == 1 :
#                  print "DEBUG_getArrayFactor (za = %.2f [deg]) : ax += %s * %s = %s" % (za[0]*(180.00/math.pi),port_current[1,i],(numpy.cos(ph[0]) + 1.0j*numpy.sin(ph[0])),ax[0])
               
        logger.info("Number of broken dipoles = %d" % count_broken )
        # set the points below the horizon to zero
        p = za >= numpy.pi/2.0
        ax[p] = 0.0
        ay[p] = 0.0
        return (ax,ay)

    def getResponse(self,az,za,freq=155e6,delays=None):
        """
        Get the full Jones matrix response of the tile including the dipole
        reponse and array factor incorporating any mutual coupling effects
        from the impedance matrix. freq in Hz.
        delays in unit steps of beamformer delays as numpy array shape (2,n_dipoles).
        az and za (radian) are numpy arrays of equal length defining a set
        of points to calculate the response for.
        Result is an array like az/za with [2][2] on the end for the Jones.
        """
        print("delays = %s , size(delays) = %d, shape(delays) = %s, n_dipoles_eda = %d" % (delays,numpy.size(delays),delays.shape,n_dipoles_eda))        
#        assert delays == None or numpy.size(delays)==2*n_dipoles_eda, "Expecting 512 delays, got %r" % str(numpy.size(delays))
        assert delays is None or (numpy.size(delays[0])==n_dipoles_eda and numpy.size(delays[1])==n_dipoles_eda), "Expecting 2 x %d delays, got %s" % ((n_dipoles_eda),delays.shape )
        print("Number of delays = %d" % numpy.size(delays))
        print("delays = %s" % (delays))
        print("n_dipoles_eda = %d" % (n_dipoles_eda))               

        # assert delays == None or numpy.size(delays)==2*n_dipoles_eda, "Expecting 512 delays, got %r" % str(numpy.size(delays))
        print("Number of delays = %d" % numpy.size(delays))
        
        if delays is None:
            delays=numpy.zeros((2,n_dipoles_eda),dtype=numpy.float32)
        print("Number of delays = %d" % numpy.size(delays))
            
        (ax,ay) = self.getArrayFactor(az,za,freq,delays)
        # get the zenith response to normalise to:
        (zax,zay) = self.getArrayFactor(numpy.array([0.0]),numpy.array([0.0]),freq) # no delays == zenith
        ax /= numpy.abs(zax)
        ay /= numpy.abs(zay)
        d = self.dipoles[0]     # for now, assume all dipoles identical FIXME
        j = d.getJones(az,za,freq)
        j[...,0,0] *= ax
        j[...,0,1] *= ax
        j[...,1,0] *= ay
        j[...,1,1] *= ay
        return j

    def getIdealDelays(self,az,za):
        """
        Calculate the ideal delays given the positions of the dipoles to point
        at the sky direction defined by az and za (in radian)
        Returns the delays for dipoles in the same order as the coords in the xpos and ypos
        """
        # calculate a unit vector to the look direction
        p_x = numpy.sin(az)*numpy.sin(za)
        p_y = numpy.cos(az)*numpy.sin(za)
        p_z = numpy.cos(za)
        unit_vec = numpy.array([p_x,p_y,p_z])
        baselines = numpy.array([self.xpos,self.ypos,self.zpos])
        print("DEBUG : getIdealDelays : using z = %s" % (self.zpos))
        delays = numpy.dot(unit_vec,baselines)
#        print "Ideal delays in nanoseconds = %s" % ((delays/vel_light)*1e9)

        delays_ns = ((delays/vel_light)*1e9)
        print("Ideal delays in nanoseconds :")
        for ant in range(0,len(delays)) :
           print("Ant_index_%03d : %.4f [ns]" % (ant,delays_ns[ant]))
           
        return delays

def convertJonesAzEl2HaDec(az,za,dec,lat):
    """
    Generate a converter for arrays of Jones matrices in Az/ZA to HA/DEC.
    This is a rotation of the orthogonal az/za unit vectors to
    the orthogonal ha/dec unit vectors by the parallactic angle
    plus some negative signs to account for the definition of
    coord axes. Multiply this by your az/za Jones to get a HA/DEC Jones.
    Inputs and outputs are numpy arrays of dimension [az][za][2][2]"""
    # parallactic angle needs HA, DEC and lat.
    (ha,dec) = h2e(az,za,lat)
    pa = calcParallacticAngle(ha,dec,lat)    
    rot = numpy.empty((ha.shape+(2,2)))
    crot = numpy.cos(pa*numpy.pi/180.0)
    srot = numpy.sin(pa*numpy.pi/180.0)
    # for clockwise rotations, use [[cos,sin],[-sin,cos]]
    rot[:,:,0,0] = crot
    rot[:,:,0,1] = srot
    rot[:,:,1,0] = -srot
    rot[:,:,1,1] = crot
    
def h2e(az,za,lat):
    """
    Horizon to equatorial.
    Convert az/za (radian) to HA/DEC (degrees, degrees)
    given an observatory latitude (degrees)
    """
    sa = numpy.sin( az )
    ca = numpy.cos( az )
    se = numpy.sin( numpy.pi/2.0 - za )
    ce = numpy.cos( numpy.pi/2.0 - za )
    sp = numpy.sin( lat*numpy.pi/180.0 )
    cp = numpy.cos( lat*numpy.pi/180.0 )

    # HA,Dec as x,y,z */
    x = - ca * ce * sp + se * cp
    y = - sa * ce
    z = ca * ce * cp + se * sp

    #To spherical */
    r = numpy.sqrt ( x*x + y*y )
    ha = numpy.arctan2 ( y, x ) * 180.0/numpy.pi
    dec = numpy.arctan2 ( z, r )* 180.0/numpy.pi
    return (ha,dec)


def e2h(ha,dec,lat):
    """
    Equatorial to horizon.
    Convert equatorial ha/dec coords (both in degs) to az,za (radian)
    given an observer latitute (degs). Returns (az,za)
    """
    ha_rad = ha*numpy.pi/180.0
    dec_rad = dec*numpy.pi/180.0
    lat_rad = lat*numpy.pi/180.0
    sh = numpy.sin(ha_rad)
    ch = numpy.cos(ha_rad)
    sd = numpy.sin(dec_rad)
    cd = numpy.cos(dec_rad)
    sp = numpy.sin(lat_rad)
    cp = numpy.cos(lat_rad)
    x = - ch * cd * sp + sd * cp;
    y = - sh * cd;
    z = ch * cd * cp + sd * sp;
    r = numpy.sqrt(x*x + y*y)
    a = numpy.atan2( y, x )
    if a < 0.0:
        az = a + 2.0*numpy.pi
    else:
        az = a
    el = numpy.atan2( z, r );
    return (az,numpi.pi/2.0 - el)


def calcParallacticAngle(ha,dec,lat):
    """
    Calculate the parallactic angle in degrees given an HA (degs)
    dec (degrees) and observatory latitude (degrees)
    """
    cl = numpy.cos(lat*numpy.pi/180.0)
    sl = numpy.sin(lat*numpy.pi/180.0)
    ch = numpy.cos(ha*numpy.pi/180.0)
    sh = numpy.sin(ha*numpy.pi/180.0)
    cd = numpy.cos(dec*numpy.pi/180.0)
    sd = numpy.sin(dec*numpy.pi/180.0)
    num = cl*sh
    den = sl*cd - cl*sd*ch
    return numpy.arctan2(num,den)*180.0/numpy.pi

def makeAZZA(npix=256):
    """
    Make azimuth and zenith angle arrays for a square image of side npix
    Projection is sine, all-sky
    Returns (az,za). Angles are in radian.
    """
    # build az and za arrays
    z = numpy.arange(npix,dtype=numpy.float32)-npix/2
    x = numpy.empty((npix,npix),dtype=numpy.float32)
    y = numpy.empty((npix,npix),dtype=numpy.float32)
    for i in range(npix):
        y[i,0:] = z
        x[0:,i] = z
    d = numpy.sqrt(x*x + y*y)/(npix/2)
    # only select pixels above horizon
    t = (d <= 1.0)
    za = numpy.ones((npix,npix),dtype=numpy.float32)*numpy.pi/2.0
    za[t] = numpy.arcsin(d[t])
    az = numpy.arctan2(y,x)
    return az,za
    
def makeUnpolInstrumentalResponse(j1,j2):
    """
    Form the visibility matrix in instrumental response from two Jones
    matrices assuming unpolarised sources (hence the brightness matrix is
    the identity matrix)
    Input: j1,j2: Jones matrices of dimension[za][az][2][2]
    Returns: [za][az][[xx,xy],[yx,yy]] where "X" and "Y" are defined by the receptors
    of the Dipole object used in the ApertureArray. Hence to get "XX", you want
    result[za][az][0][0] and for "YY" you want result[za][az][1][1]
    """
    result = numpy.empty_like(j1)

    result[:,:,0,0] = j1[:,:,0,0]*j2[:,:,0,0].conjugate() + j1[:,:,0,1]*j2[:,:,0,1].conjugate()
    result[:,:,1,1] = j1[:,:,1,0]*j2[:,:,1,0].conjugate() + j1[:,:,1,1]*j2[:,:,1,1].conjugate()
    result[:,:,0,1] = j1[:,:,0,0]*j2[:,:,1,0].conjugate() + j1[:,:,0,1]*j2[:,:,1,1].conjugate()
    result[:,:,1,0] = j1[:,:,1,0]*j2[:,:,0,0].conjugate() + j1[:,:,1,1]*j2[:,:,0,1].conjugate()
    return result

def makePolInstrumentResponse(j1,j2,b):
    """
    Form the instrument response from two Jones matrices with an
    arbitrary source brightness matrix, hence arbitrary polarisation
    Returns: (xx,yy,xy,yx) where "X" and "Y" are defined by the receptors
    of the Dipole object used in the ApertureArray
    """
    # FIXME: need to work out how to do this in vectorised way.

def plotDipoleJones(d,freq=155e6):
    """
    Utility to make a plot of a dipole Jones matrix for debugging
    """

    import matplotlib.pyplot as plt
    (az,za) = makeAZZA()
    logger.info("plotting dipole Jones response for type: "+d.type+', freq (MHz): '+str(freq/1e6))
    j = d.getJones(az,za,freq)

    plt.imshow(j[:,:,0,0].real)
    plt.title('EDA '+str(freq/1e6)+'MHz dipole J00 voltage real')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J00_voltage_re_'+str(freq/1e6)+'MHz_dipole.png')
    plt.clf()

    plt.imshow(j[:,:,0,1].real)
    plt.title('EDA '+str(freq/1e6)+'MHz dipole J01 voltage real')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J01_voltage_re_'+str(freq/1e6)+'MHz_dipole.png')
    plt.clf()

    plt.imshow(j[:,:,1,0].real)
    plt.title('EDA '+str(freq/1e6)+'MHz dipole J10 voltage real')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J10_voltage_re_'+str(freq/1e6)+'MHz_dipole.png')
    plt.clf()

    plt.imshow(j[:,:,1,1].real)
    plt.title('EDA '+str(freq/1e6)+'MHz dipole J11 voltage real')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J11_voltage_re_'+str(freq/1e6)+'MHz_dipole.png')
    plt.clf()


def plotArrayJones(j,freq,za):
    """
    Utility to plot the output of tile Jones matrices
    """
    import matplotlib.pyplot as plt

    plt.imshow(numpy.abs(j[:,:,0,0]))
    plt.title('EDA '+str(freq/1e6)+'MHz J00 voltage mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J00_voltage_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.abs(j[:,:,0,1]))
    plt.title('EDA '+str(freq/1e6)+'MHz J01 voltage mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J01_voltage_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.abs(j[:,:,1,0]))
    plt.title('EDA '+str(freq/1e6)+'MHz J10 voltage mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J10_voltage_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.abs(j[:,:,1,1]))
    plt.title('EDA '+str(freq/1e6)+'MHz J11 voltage mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_J11_voltage_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

def plotVisResponse(j,freq,za):
    """
    Utility to plot the visibilty XX,YY,XY and YX response of the array for
    an unpolarised 1Jy source
    Input: j a visibility matrix (complex) of dimensions [za][az][2][2]
    """
    import matplotlib.pyplot as plt

    vis = makeUnpolInstrumentalResponse(j,j)
    plt.imshow(numpy.abs(vis[:,:,0,0]))
    plt.title('EDA '+str(freq/1e6)+'MHz XX mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_XX_mag_'+str(freq/1e6)+'MHz_ZA='+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.abs(vis[:,:,1,1]))
    plt.title('EDA '+str(freq/1e6)+'MHz YY mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_YY_mag_'+str(freq/1e6)+'MHz_ZA='+str(za)+'.png')
    plt.clf()
    
    plt.imshow(numpy.abs(vis[:,:,0,1]))
    plt.title('EDA '+str(freq/1e6)+'MHz XY mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_XY_mag_'+str(freq/1e6)+'MHz_ZA='+str(za)+'.png')
    plt.clf()
    
    plt.imshow(numpy.abs(vis[:,:,1,0]))
    plt.title('EDA '+str(freq/1e6)+'MHz YX mag ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_YX_mag_'+str(freq/1e6)+'MHz_ZA='+str(za)+'.png')
    plt.clf()
    

def plotArrayFactors(ax,ay,za):
    """
    Utility to make a plot of an array factor for debugging
    """
    import matplotlib.pyplot as plt

    plt.imshow(numpy.abs(ax))
    plt.title('EDA '+str(freq/1e6)+'MHz X voltage mag array factor. ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_array_voltage_X_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.angle(ax))
    plt.title('EDA '+str(freq/1e6)+'MHz X voltage phase array factor ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_array_voltage_X_ph_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.abs(ay))
    plt.title('EDA '+str(freq/1e6)+'MHz Y voltage mag array factor ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_array_voltage_Y_mag_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

    plt.imshow(numpy.angle(ay))
    plt.title('EDA '+str(freq/1e6)+'MHz Y voltage phase array factor ZA='+str(za))
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.savefig('/tmp/EDA_array_voltage_Y_ph_'+str(freq/1e6)+'MHz_ZA'+str(za)+'.png')
    plt.clf()

def addWCStoFITS(filename,freq_MHz,lst_hours):
    dec=-26.7033    # MWA latitude
    point_az = 0.0
    point_za = 0.0
    hdulist = pyfits.open(filename, mode='update')
    prihdr = hdulist[0].header
    siz = prihdr['NAXIS1']
    pixscale=180.0/siz # for all-sky
    prihdr.set('CRPIX1',siz/2+1)
    prihdr.set('CDELT1',pixscale/numpy.pi)
    prihdr.set('CRVAL1',lst_hours*15.0)
    prihdr.set('CTYPE1',"RA---SIN")
    prihdr.set('CRPIX2',siz/2+1)
    prihdr.set('CDELT2',pixscale/numpy.pi)
    prihdr.set('CRVAL2',dec)
    prihdr.set('CTYPE2',"DEC--SIN")
    prihdr.set('BEAM_AZ',point_az,' [deg]')
    prihdr.set('BEAM_ZA',point_za,' [deg]')
    prihdr.set('FREQ',freq_MHz*1e6,' [Hz]')
    hdulist.close()

def save_feko_file( za_map, az_map, vis, freq, out_filename="feko.txt" ) :
    assert az_map.shape==za_map.shape, "Input az and za arrays must have same dimenions"    
    out_file=open(out_filename,"w")

    # HEADER:
    outline="#TEST \n"
    out_file.write(outline)
    
#    shape = az_map.
    size_x = az_map.shape[0]
    size_y = az_map.shape[1]
    beam_power = numpy.abs(vis[:,:,0,0]).copy()
    beam_power[numpy.isnan(beam_power)]=0.00    
    za_map[numpy.isnan(za_map)]=-1000
    for y in range(0,size_y):
       for x in range(0,size_x):
          za=za_map[x,y]*(180.00/numpy.pi)
          az=az_map[x,y]*(180.00/numpy.pi)
          
          if za >= 0 :
             outline = "%.8f\t%.8f\t0\t0\t0\t0\t0\t0\t%.20f\n" % (za,az,beam_power[x,y])
             out_file.write(outline)

    out_file.close() 
    
    za_file_name="za.fits"
    pyfits.writeto(za_file_name,za_map,clobber=True)        

    az_file_name="az.fits"
    pyfits.writeto(az_file_name,az_map,clobber=True)        
    print("saved ???")


def exportToFITS(j,freq,za,lst_hours,za_map=None,az_map=None,save_feko=False):
    vis = makeUnpolInstrumentalResponse(j,j)
    filename = "EDA_beam_ZA_"+str(za)+"_"+str(freq/1e6)+"MHz_XX.fits"
    try:
        os.remove(filename)
    except:
        sys.exc_clear()
    pyfits.writeto(filename, numpy.abs(vis[:,:,0,0]))
    addWCStoFITS(filename,freq/1e6,lst_hours)
    filename = "EDA_beam_ZA_"+str(za)+"_"+str(freq/1e6)+"MHz_YY.fits"
    try:
        os.remove(filename)
    except:
        sys.exc_clear()
    pyfits.writeto(filename, numpy.abs(vis[:,:,1,1]))
    addWCStoFITS(filename,freq/1e6,0.0)
    
    if save_feko and za_map is not None and az_map is not None :
       save_feko_file( za_map, az_map, vis, freq )          


tiles_initialised=False
d=None
d_ideal=None
dipoles = []
dipoles_ideal = []
tile=None
tile_ideal=None
g_freq=-100

def init_tiles( gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00, doplots=False, freq=240e6, type='short', xpos=None, ypos=None, zpos=None, station_name="EDA", projection="zea" ) :
    global tiles_initialised
    global d
    global d_ideal
    global dipoles
    global dipoles_ideal
    global tile
    global tile_ideal
    global g_freq
    
    if freq != g_freq :
       print("init_tiles(%.2f,%.2f): frequency change %.2f -> %.2f" % (gain_sigma_dB,gain_sigma_ph_160mhz,g_freq,freq))
       logger.debug("init_tiles(%.2f,%.2f): frequency change %.2f -> %.2f" % (gain_sigma_dB,gain_sigma_ph_160mhz,g_freq,freq))
       tiles_initialised = False

    if not tiles_initialised :
       logger.debug("init_tiles(%.2f,%.2f): initialising objects" % (gain_sigma_dB,gain_sigma_ph_160mhz))
       print("init_tiles(%.2f,%.2f): initialising objects" % (gain_sigma_dB,gain_sigma_ph_160mhz))
       d = Dipole(type=type,station_name=station_name,projection=projection)
       d_ideal = Dipole(type=type,station_name=station_name,projection=projection)
       if doplots:
           plotDipoleJones(d,freq)
       dipoles = []
       for i in range(n_dipoles_eda):
           dipoles.append(Dipole(type,station_name=station_name))

       dipoles_ideal = []
       for i in range(n_dipoles_eda):
           dipoles_ideal.append(d_ideal)

       # randomise the gains on the dipoles. MWA dipoles have gain amplitude of 19dB with routhly 1 sigma variation about 1 dB.
       # Gain measurements :
       # Dave Kenney   : 0.5 or 0.7 for beamformes
       # Abraham Neben : ~0.4dB beamformers 
       # Short EDA-dipole test ( eda_short_dipole_tests_201611.odt ) : ~0.5 dB 
       # Seems that good value is between 0.5-1.0 dB variation
       # Brian's boards measurements ?
       gain_sigma_ph = gain_sigma_ph_160mhz * freq/160e6 # based on DaveK's measurements
       for dipole in dipoles:
           newgain = numpy.power(10,gain_sigma_dB*numpy.random.randn(2)/10)
           phase = gain_sigma_ph * numpy.random.randn(1) * numpy.pi/180.0
           phase_c = numpy.complex(numpy.cos(phase),numpy.sin(phase))
           dipole.gain[0,0] *= newgain[0]*phase_c
           dipole.gain[1,1] *= newgain[1]*phase_c
       tile = ApertureArray( dipoles=dipoles, xpos=xpos, ypos=ypos, zpos=zpos )
       tile_ideal = ApertureArray( dipoles=dipoles_ideal, xpos=xpos, ypos=ypos, zpos=zpos )
  
       tiles_initialised = True
       g_freq = freq
    else :
       print("init_tiles(%.2f,%.2f) : Already called - just returning old objects" % (gain_sigma_dB,gain_sigma_ph_160mhz))
       logger.debug("init_tiles(%.2f,%.2f) : Already called - just returning old objects" % (gain_sigma_dB,gain_sigma_ph_160mhz))
       
    return (d,d_ideal,dipoles,dipoles_ideal,tile,tile_ideal)       


# execute a series of tests if invoked from the command line
# was pointing_za_deg=6.81 
def get_single_dipole_beam( za, az, pointing_za_deg=0.00, pointing_az_deg=0.00, resolution=512, delays=None, zenithnorm=True, power=True, jones=False, freq = 240e6, 
                            lst=0.00, gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00, dipole_type='short', 
                            xpos=None, ypos=None, zpos=None,  # list of antenna positions can overwrite the default list 
                            use_beam_fits=True, station_name="EDA", projection="zea"
                          ) :    
    if use_beam_fits :
       dipole_type = "fits_beam"
    else :
       logger.warning("get_single_dipole_beam : only implemented to use beam FITS files")    
    logger.setLevel(logging.DEBUG)
    logger.info("get_single_dipole_beam , use_beam_fits=%s, dipole_type = %s, projection = %s" % (use_beam_fits,dipole_type,projection))
    
    az_deg = az*(180.00/math.pi)
    za_deg = za*(180.00/math.pi)
    freq_mhz = freq/1e6
    
    (beam_x,current_fits_beam_x,x_pixel_x,y_pixel_x) = fits_beam.get_fits_beam_multi( az_deg , za_deg , freq_mhz, polarisation='X', station_name=station_name, power=True, projection=projection )  # get_fits_beam
    (beam_y,current_fits_beam_y,x_pixel_y,y_pixel_y) = fits_beam.get_fits_beam_multi( az_deg , za_deg , freq_mhz, polarisation='Y', station_name=station_name, power=True, projection=projection )  # get_fits_beam
    
    return (beam_x,beam_y)
    


# execute a series of tests if invoked from the command line
# was pointing_za_deg=6.81 
def get_eda_beam( za, az, pointing_za_deg=0.00, pointing_az_deg=0.00, resolution=512, delays=None, zenithnorm=True, power=True, jones=False, freq = 240e6, 
                  lst=0.00, gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00, dipole_type='short', 
                  xpos=None, ypos=None, zpos=None,  # list of antenna positions can overwrite the default list 
                  use_beam_fits=False, station_name="EDA", projection="zea"
                ) :    
    if use_beam_fits :
       dipole_type = "fits_beam"
    logger.setLevel(logging.DEBUG)
    logger.info("get_eda_beam , use_beam_fits=%s, dipole_type = %s" % (use_beam_fits,dipole_type))
    
    doplots=False
    exportBeam=True
    lat = -26.7
#    (az,za) = makeAZZA(npix=resolution)

    # parallactic angle
    (ha,dec) = h2e(az,za,lat)
    pa = calcParallacticAngle(ha,dec,lat)

    (d,d_ideal,dipoles,dipoles_ideal,tile,tile_ideal) = init_tiles( gain_sigma_dB=gain_sigma_dB , gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, doplots=doplots, freq=freq, type=dipole_type, xpos=xpos, ypos=ypos, zpos=zpos, station_name=station_name, projection=projection  )
#    d = Dipole()
#    d_ideal = Dipole()
#    if doplots:
#        plotDipoleJones(d,freq)
#    dipoles = []
#    for i in range(n_dipoles_eda):
#        dipoles.append(Dipole('short'))
#
#    dipoles_ideal = []
#    for i in range(n_dipoles_eda):
#        dipoles_ideal.append(d_ideal)
#
#    # randomise the gains on the dipoles. MWA dipoles have gain amplitude of 19dB with routhly 1 sigma variation about 1 dB.
#    gain_sigma_ph = gain_sigma_ph_160mhz * freq/160e6 # based on DaveK's measurements
#    for dipole in dipoles:
#        newgain = numpy.power(10,gain_sigma_dB*numpy.random.randn(2)/10)
#        phase = gain_sigma_ph * numpy.random.randn(1) * numpy.pi/180.0
#        phase_c = numpy.complex(numpy.cos(phase),numpy.sin(phase))
#        dipole.gain[0,0] *= newgain[0]*phase_c
#        dipole.gain[1,1] *= newgain[1]*phase_c
#    tile = ApertureArray(dipoles=dipoles)
#    tile_ideal = ApertureArray(dipoles=dipoles_ideal)

    print("pointing_az_deg = %.8f , pointing_za_deg = %.8f" % (pointing_az_deg,pointing_za_deg))
    delays_0_75_ideal = tile.getIdealDelays( float(pointing_az_deg)*(numpy.pi/180.00) , float(pointing_za_deg)*(numpy.pi/180.00) )/DQ

    # systematic error in delays :
    # delays_0_75_ideal = delays_0_75_ideal*1.05
    
#    delays_0_75_ideal = tile.getIdealDelays(0,numpy.pi*float(6.81)/float(180))/DQ
    delays_0_75_quantised = numpy.round(delays_0_75_ideal)
    delays1=numpy.array([[0]*n_dipoles_eda,[0]*n_dipoles_eda],dtype=numpy.float32)
    za_delays = {'0':delays1*0,'15':numpy.array([delays_0_75_ideal,delays_0_75_ideal])}

    (ax0,ay0) = tile.getArrayFactor(az,za,freq,za_delays['15'])
    val=numpy.abs(ax0)
    val_no_nan=val.copy()
    numpy.isfinite(val, val_no_nan)    
#    val_finite=val[val_no_nan]
#    val_max=numpy.amax(val_finite)
    val_max=numpy.nanmax(val)
    if val.shape[0]>1 and val.shape[1]>1 :
       print("VALUE : %.8f %.8f %.8f" % (freq,val_max,val[resolution/2,resolution/2]))                

    if doplots:
        for za_delay in za_delays:
            logger.debug("ZA is: %s. Delays are: %r" % (za_delay,za_delays[za_delay]))
            (ax0,ay0) = tile.getArrayFactor(az,za,freq,za_delays[za_delay])
            logger.info("plotting Array factor voltage for ZA "+str(za_delay))
            plotArrayFactors(ax0,ay0,za_delay)
            logger.info("plotting tile Jones response for ZA "+str(za_delay))
            j = tile_ideal.getResponse(az,za,freq,za_delays[za_delay])
            plotArrayJones(j,freq,za_delay)
            logger.info("Plotting visbility response for two identical tiles ZA "+str(za_delay))
            plotVisResponse(j,freq,za_delay)

    myza='15' # set to 0 for zenith  '15' to use pointing as specified by parameters
    za_delay = za_delays[myza]
    j = tile.getResponse(az,za,freq,za_delay)
    if exportBeam:
        # export the beam model
        exportToFITS(j,freq,myza,lst,za,az)

        # export the sky area weights for a sine projection image
        (k,l) = numpy.where(za < 0.99*numpy.pi/2)
        area_scale_factor = za*0.0
        area_scale_factor[k,l] = 1.0/numpy.cos(za[k,l])
        fname='sky_scale.fits'
        try:
            os.remove(fname)
        except:
            sys.exc_clear()
        pyfits.writeto(fname, area_scale_factor)

    vis = makeUnpolInstrumentalResponse(j,j)
    if not power:
       return (numpy.sqrt(vis[:,:,0,0].real),numpy.sqrt(vis[:,:,1,1].real))
    else:
       return (vis[:,:,0,0].real,vis[:,:,1,1].real)

def shift_antenna( ant_idx, dx=0, dy=0 ) :
   global xs
   global yx
   
   xs[ant_idx] =  xs[ant_idx] + dx
   ys[ant_idx] =  ys[ant_idx] + dy


if __name__ == "__main__":
   logger.setLevel(logging.DEBUG)
   doplots=False
   exportBeam=True
   lat = -26.7
   freq = 240e6
   
   pointing_az_deg=0
   if len(sys.argv) >= 2:
      pointing_az_deg=float(sys.argv[1])

   pointing_za_deg=0
   if len(sys.argv) >= 3:
      pointing_za_deg=float(sys.argv[2])
      

   freq=100e6
   if len(sys.argv) >= 4:
      freq=float(sys.argv[3])

   gain_sigma_dB=0.7
   gain_sigma_ph_160mhz=3
      
   npix=64 # was : 512   
   (az,za) = makeAZZA(npix=npix)
   out_f_za = open("out_za.txt", "w" )
   out_f_az = open("out_az.txt", "w" )
   
   for y in range(npix-1,-1,-1):
      line_az = ""
      line_za = ""
      for x in range(0,npix) :
         az_str = "%06.1f " % (az[x,y]*(180.00/math.pi))
         za_str = "%06.1f " % (za[x,y]*(180.00/math.pi))
         line_az += az_str
         line_za += za_str
         
      out_f_az.write( line_az + "\n" )
      out_f_za.write( line_za + "\n" )               
           
   out_f_za.close()
   out_f_az.close()



   # parallactic angle
   (ha,dec) = h2e(az,za,lat)
   pa = calcParallacticAngle(ha,dec,lat)
                                                                        
   (beam_x,beam_y) = get_eda_beam( za, az, pointing_az_deg=pointing_az_deg, pointing_za_deg=pointing_za_deg, freq=freq, gain_sigma_dB=gain_sigma_dB, gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, resolution=npix, use_beam_fits=False )     
#   print("beam_x = %.8f" % (beam_x[0,0]))
#   print("beam_y = %.8f" % (beam_y[0,0]))
                                                           