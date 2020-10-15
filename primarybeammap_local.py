#!/usr/bin/env python
"""
primarybeammap.py --freq=202.24 --beamformer=0,0,0,1,3,3,3,3,6,6,6,6,8,9,9,9 --gps=20110926210616

main task is:
make_primarybeammap()


"""
from __future__ import print_function
# from mwapy import ephem_utils
# import primary_beam 
import beam_tools
import eda_beam
import read_feko
import fits_beam

import sys
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy,math
import os
import logging
import matplotlib
if not 'matplotlib.backends' in sys.modules:
    matplotlib.use('agg')
import matplotlib.pyplot as pylab
# import pylab
from scipy.interpolate import RegularGridInterpolator

import astropy
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import numpy as np

EPS=numpy.finfo(numpy.float64).eps #machine epsilon

defaultcolor='k'
defaultsize=8
contourlevels=[0.01, 0.1, 0.25, 0.5, 0.75]
file_haslam = None

# configure the logging
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('primarybeammap')
logger.setLevel(logging.WARNING)

# radio_image='radio408.RaDec.fits'

MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)   

add_sun=False
sun_ra=None  # test :290.46450056
sun_dec=None # test : -26.70331900
# SUN :
#   double sun_T_b = 600000.00; // from Claude Mercier and Gilbert Chambe 2009 paper 
#   double sun_radius_arcsec = 1818.00/2.00;
# Lantos & Avignon - Sun ~ 60000Jy 
sun_T_b = 1000000.00 # Lantos et all 1975 etc - difficult to dind one good value 
sun_radius_arcsec = (31.00/2.00)*60.00 # Diameter ~ 31 arcmin 
sun_radius_deg = sun_radius_arcsec / 3600.00



def set_sun( add_sun_param, sun_ra_param, sun_dec_param, scale_size=1.00 ) :
   global add_sun
   global sun_ra
   global sun_dec
   global sun_radius_arcsec
   global sun_radius_deg   
 
   add_sun = add_sun_param
   sun_ra  = sun_ra_param
   sun_dec = sun_dec_param
   
   sun_radius_arcsec = sun_radius_arcsec * scale_size
   sun_radius_deg = sun_radius_arcsec / 3600.00
   
   print("Sun parameters radius = %.2f [arcsec] and coordinates will be calculated for each timestamp (new version)" % (sun_radius_arcsec))

def shift_antenna( ant_idx, dx=0, dy=0 ) :
   eda_beam.shift_antenna( ant_idx=ant_idx, dx=dx, dy=dy )   
   

######################################################################
def get_azza_arrays_fov(gridsize=361,fov=180.0):
    """
    Converted from Randall Wayth's IDL code. 
    
    Make Az,ZA arrays for a field of view (degrees)
    az array value range is -180 to 180 degrees

    gridsize is the grid size along an edge
    fov=180 degrees is the visible field of view
    """
    if fov > 180.0:
        logger.error("FOV of %s is too large. Max: 180 deg" % (fov))
        return None
    mask = numpy.zeros((gridsize,gridsize),float)
    za_grid = numpy.zeros((gridsize,gridsize),float)
    #create u,v plane (given as c, r)
    a=numpy.arange(-1,1+1.0/(gridsize-1),2.0/(gridsize-1))
    c,r=numpy.meshgrid(a,a)
    myfov = math.sin(fov/2*math.pi/180)
    dsqu = (c**2+r**2)*(myfov)**2
    p = (dsqu<(1.0+EPS))
    za_grid[p]=numpy.arcsin(dsqu[p]**0.5)
    print('Using standard orthographic projection')
    az_grid=numpy.arctan2(c,r)
    mask[p] = 1.0 #set mask 
    p = dsqu>=(1.0+EPS)
    za_grid[p] = math.pi/2.0 #set ZA outside of fov to 90 deg
    return az_grid*180.0/math.pi, za_grid*180.0/math.pi

def ang_dist( ra1_deg, dec1_deg, ra2_deg, dec2_deg ):
   ra1=ra1_deg*(math.pi/180.00)
   dec1=dec1_deg*(math.pi/180.00)
   ra2=ra2_deg*(math.pi/180.00)
   dec2=dec2_deg*(math.pi/180.00)
   
   cos_value = math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos( ra1 - ra2 )
   dist_value = math.acos( cos_value )
   return dist_value
                

######################################################################
#def map_sky(skymap,lst,lat,az_grid,za_grid):
#    """
#    Converted from Randall Wayth's IDL code.
#
#    Map skymap onto grid of arbitrary size"""
#    out=az_grid*0.0 # new array for gridded sky

#    ha_grid,dec_grid=ephem_utils.horz2eq(az_grid,90-za_grid,lat) # get grid in ha, dec 
#    size_dec=skymap.shape[0]
#    size_ra=skymap.shape[1]
#    p = za_grid<90.0+EPS #array indices for visible sky

    # the following assumes RA=0 in centre
    # of the sky image and increases to the left.
#    ra = (lst - ha_grid/15.0) %  24.0
#    ra_index = (((36-ra) % 24)/24)*size_ra
#    dec_index = (dec_grid/180.0+0.5)*size_dec

#    print ra_index.min(), ra_index.max()
#    print dec_index.min(), dec_index.max()

    #select pixels of sky map, using ra and dec index values
    # rounded down to nearest index integer
    #print p
    #print numpy.rint(ra_index[p]),numpy.rint(dec_index[p])
#    print numpy.rint(ra_index[p]).astype(int)
#    print numpy.rint(dec_index[p]).astype(int)
#    print skymap.shape
#    out[p]=skymap[dec_index[p].astype(int),ra_index[p].astype(int)] 
#    return out

def map_sky_astropy(skymap,RA,dec,gps,az_grid,za_grid,epoch='B1950'):
    """Reprojects Haslam map onto an input az, ZA grid.
Inputs:
skymap 
RA - 1D range of RAs (deg)
dec - 1D range of decs (deg)
gps - GPS time of observation
az_grid - grid of azes onto which we map sky
za_grid - grid of ZAs onto which we map sky

    """
    #Get az, ZA grid transformed to equatorial coords
    grid2eq=horz2eq(az_grid, za_grid, gps, epoch=epoch)
    print('grid2eq',grid2eq['RA'].shape)
   

    #Set up interp function using sky map
    #flip so interpolation has increasing values
    #TODO: I don't think this will affect outcome! 
    my_interp_fn=RegularGridInterpolator((dec, RA[::-1]), skymap[:,::-1],fill_value=None)
    #fill_value = None means that values outside domain are extrapolated.
    #fill_value=nan would be preferable, but this causes error due to bug in scipy<15.0, as per
    #https://github.com/scipy/scipy/issues/3703
    
    #interpolate map onto az,ZA grid
    print(np.min(grid2eq['dec']), np.max(grid2eq['dec']))
    print(np.min(grid2eq['RA']), np.max(grid2eq['RA']))
    #Convert to RA=-180 - 180 format (same as Haslam)
    #We do it this way so RA values are always increasing for RegularGridInterpolator
    grid2eq['RA'][grid2eq['RA']>180]=grid2eq['RA'][grid2eq['RA']>180]-360
    
    print(np.min(grid2eq['dec']), np.max(grid2eq['dec']))
    print(np.min(grid2eq['RA']), np.max(grid2eq['RA']))
    my_map=my_interp_fn(np.dstack([grid2eq['dec'], grid2eq['RA']]))
#    print "np.vstack([grid2eq['dec'], grid2eq['RA']])",np.vstack([grid2eq['dec'], grid2eq['RA']]).shape
#    print "np.hstack([grid2eq['dec'], grid2eq['RA']])",np.hstack([grid2eq['dec'], grid2eq['RA']]).shape
#    print "np.dstack([grid2eq['dec'], grid2eq['RA']])",np.dstack([grid2eq['dec'], grid2eq['RA']]).shape    

    # test adding sun
    if add_sun : 
       # 2020-10-13 - always calculate current Sun position :
       # uxtime = Time( gps, scale='gps', format="unix" )       
       uxtime = gps + 315964783 # TODO : should use proper function, but cannot find it now ...
       (sun_ra,sun_dec, sun_az, sun_alt, sun_za ) = fits_beam.sun_position( unix_time=uxtime )
    
       if sun_alt > 0 :
          print("Adding sun at (ra,dec)=(%.4f,%.4f) [deg] (az,elev)=(%.4f,%.4f) [deg] for gps=%d / uxtime = %d" % (sun_ra,sun_dec,sun_az,sun_alt,gps,uxtime))
          print("dec shape = %d x %d" % (grid2eq['dec'].shape[0],grid2eq['dec'].shape[1]))
          for dy in range(0,grid2eq['dec'].shape[1]):
             for dx in range(0,grid2eq['dec'].shape[0]):
                dec_px=grid2eq['dec'][dx,dy]
                ra_px=grid2eq['RA'][dx,dy]
          
                dist_arcsec = ang_dist(ra_px,dec_px,sun_ra,sun_dec)*(180.00/math.pi)*3600.00
                if dist_arcsec < sun_radius_arcsec :
                   my_map[dx,dy] = sun_T_b
                   print("Sun detected")
       else : 
          print("WARNING : Sun below horizon at (ra,dec)=(%.4f,%.4f) [deg] (az,elev)=(%.4f,%.4f) [deg] for gps=%d / uxtime = %d" % (sun_ra,sun_dec,sun_az,sun_alt,gps,uxtime))
       
        
    return my_map
     
def eq2horz(ra, dec, gps): 
    """Convert from equatorial (RA, dec) to horizontal (az, ZA)"
Returns Az (CW from North) and ZA in degrees at a given time,
Inputs: 
time - GPS time"""
    
    coords = SkyCoord(ra=ra, dec=dec,equinox='J2000',unit=astropy.units.deg)    
    coords.location=MWA_POS

    #convert GPS to an astropy time object        
    coords.obstime=Time(gps,format='gps',scale='utc')
    #get sidereal_time to reduced precision, by explicitly setting the offset of UT1 from UTC:
    #or the current time you can also get it much more precisely following the instructions in http://docs.astropy.org/en/latest/time/index.html#transformation-offsets)
    logger.warning('Using approximate sidereal time:') 
    coords.obstime.delta_ut1_utc = 0

    
    logger.info('Calculating az, ZA at time %s', coords.obstime) 
    mycoords=coords.transform_to('altaz')
    return {'Az':mycoords.az.deg,'ZA':90-mycoords.alt.deg}

g_printed_info = False

# J2000 = ICRS frame
# B1950 = FK4 frame
def horz2eq(az, ZA, gps, epoch='B1950' ):  #  does not seem to matter if B1950
    """Convert from horizontal (az, ZA) to equatorial (RA, dec)"
Returns RA, dec,
Inputs: 
time - GPS time"""
    global g_printed_info

    if not g_printed_info :
       print("horz2eq ( epoch = %s )" % (epoch))
       g_printed_info = True


    time=Time(gps,format='gps',scale='utc')
#    logger.info('Calculating az, ZA at time %s', coords.obstime) 
    coords=SkyCoord(alt = 90-ZA, az = az, obstime = time,
                    frame = 'altaz',unit=astropy.units.deg,equinox=epoch,
                    location=MWA_POS)

    #convert GPS to an astropy time object        
#    coords.obstime=Time(gps,format='gps',scale='utc')

    # MS : warning 20191025 B1950 (HASLAM MAP) = FK4 :
    if epoch == 'B1950' :
       return {'RA':coords.fk4.ra.deg,'dec':coords.fk4.dec.deg}
    
    if epoch == 'J2000' :
       return {'RA':coords.icrs.ra.deg,'dec':coords.icrs.dec.deg}


    print("ERROR : unknown epoch = %s" % (epoch))
    
    return {None,None}

def get_Haslam(freq, dir=None, scaling=-2.55,radio_image='radio408.RaDec.fits'):
    """get the Haslam 408 MHz map.
Outputs
RA - RA in degrees (-180 - 180)
dec - dec in degrees""" 
    global file_haslam 

    dir=os.path.dirname(__file__)
    if (len(dir)==0):
        dir='.'
    radio_image_touse=dir + '/' + radio_image
    
    if file_haslam is None :
       logger.info("Loading 408 MHz map from %s..." % radio_image_touse)
       print("Loading 408 MHz map from %s..." % radio_image_touse)
       #radio_image_touse='/data/das4/packages/MWA_Tools/mwapy/pb/radio408.RaDec.fits' 
       #radio_image_touse=radio_image
       if not os.path.exists(radio_image_touse):
           logger.error("Could not find 408 MHz image: %s\n" % (radio_image_touse))
           return None
       try:
           logger.info("Loading 408 MHz map from %s..." % radio_image_touse)
           file_haslam = pyfits.open(radio_image_touse)
       except:
           logger.error("Error opening 408 MHz image: %s\n" % (radio_image_touse))
           return None
    else : 
       logger.info("File %s already opened (re-using)\n" % (radio_image_touse))
       print("File %s already opened (re-using)" % (radio_image_touse)) 
        
    skymap=file_haslam[0].data[0]/10.0 #Haslam map is in 10xK
    skymap=skymap*(freq/408.0e6)**scaling #Scale to frequency
    
    RA_1D=(file_haslam[0].header.get('CRVAL1')+(numpy.arange(1,skymap.shape[1]+1)-file_haslam[0].header.get('CRPIX1'))*file_haslam[0].header.get('CDELT1'))#/15.0
    dec_1D=file_haslam[0].header.get('CRVAL2')+(numpy.arange(1,skymap.shape[0]+1)-file_haslam[0].header.get('CRPIX2'))*file_haslam[0].header.get('CDELT2')

    return {'skymap':skymap, 'RA':RA_1D,'dec':dec_1D} #RA, dec in degs
    
           ######################################################################
           
def get_LST(gps):
    time=Time(gps,format='gps',scale='utc')
    time.delta_ut1_utc = 0.
    LST=time.sidereal_time('apparent',MWA_POS.longitude.value)
    return LST.value #keep as decimal hr
    


#FIXME: most the arguments in make_primarybeammap are not needed
def make_primarybeammap(gps, delays, frequency, model, beams=None, extension='png', 
                        plottype='beamsky',figsize=14,title=None, directory=None,resolution=1000,out_filename="out.txt",
                        pointing_za_deg=0.00, pointing_az_deg=0.00,gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00, dipole_type='short',
                        radio_image='radio408.RaDec.fits', T_rcv=0, rcv_noise_coupling=0.00, use_trcv=False, xpos=None, ypos=None, zpos=None,
                        feko_file=None, interpolate=False, interpolate_to_nearest=False, pol_list_string="XX",
                        use_beam_fits=False, station_name="EDA", mean_dipole_beam=False,
                        projection="zea" ): 
    """
    """                
    print("Output beam file resolution = %d , output directory = %s , projection = %s" % (resolution,directory,projection))
#    (az_grid, za_grid) = beam_tools.makeAZZA(resolution,'ZEA') #Get grids in radians
    za_grid = None
    az_grid = None
    n_total = 0 
    dOMEGA  = 0.00
    
    if use_beam_fits :
        # !!! may require using ZEA - is it matters at all ...
        (az_grid, za_grid,n_total,dOMEGA) = fits_beam.makeAZZA(resolution,projection='ZEA') # 'ZEA' works fine, but is incompatible with the beam projection - check with RW 
    else :    
        (az_grid, za_grid,n_total,dOMEGA) = beam_tools.makeAZZA_dOMEGA(resolution,'ZEA') # this one is totally different 
    az_grid=az_grid*180/math.pi
    za_grid=za_grid*180/math.pi
    #az_grid+=180.0
    alt_grid=90-(za_grid)
    lst=get_LST(gps)

    # first go from altitude to zenith angle
    theta=(90-alt_grid)*math.pi/180
    phi=az_grid*math.pi/180

    # beams={}
    if beams is None : 
        beams={}

        if feko_file is not None : 
            print("INFO : beams = None and feko_file = %s -> reading beam from FEKO file -> may take a while ..." % (feko_file))
            beams['XX'] = read_feko.get_feko_beam( theta, phi, feko_file=feko_file, interpolate=interpolate, interpolate_to_nearest=interpolate_to_nearest )
            beams['YY'] = beams['XX'].copy()
        else :
            if mean_dipole_beam :
               print("INFO : calculating single dipole %s beam , projection = %s" % (station_name,projection))
               beams['XX'],beams['YY']=eda_beam.get_single_dipole_beam( theta, phi, resolution=resolution, zenithnorm=True, power=True, freq=frequency, lst=lst, 
                                                                        pointing_za_deg=pointing_za_deg, pointing_az_deg=pointing_az_deg, gain_sigma_dB=gain_sigma_dB, 
                                                                        gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, dipole_type=dipole_type, xpos=xpos, ypos=ypos, zpos=zpos,
                                                                        use_beam_fits=use_beam_fits, station_name=station_name, projection=projection )
            else :
               print("INFO : beams = None and feko_file = None -> generating normal EDA beam ( projection = %s )" % (projection))
               beams['XX'],beams['YY']=eda_beam.get_eda_beam( theta, phi, resolution=resolution, zenithnorm=True, power=True, freq=frequency, lst=lst, 
                                                              pointing_za_deg=pointing_za_deg, pointing_az_deg=pointing_az_deg, gain_sigma_dB=gain_sigma_dB, 
                                                              gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, dipole_type=dipole_type, xpos=xpos, ypos=ypos, zpos=zpos,
                                                              use_beam_fits=use_beam_fits, station_name=station_name, projection=projection )
    else :
        print("INFO : beams != None -> using externally provided beam")
                                                                 
#    pols=['XX','YY']
    pols=pol_list_string.split(",") # pols=['XX']
    
    # find beam maximum :
    max_beam_value = -1000
    max_beam_value_xx = -1
    max_beam_value_yy = -1
    for xx in range(beams['XX'].shape[0]) :
       for yy in range(beams['XX'].shape[1]) :
          if beams['XX'][xx,yy] > max_beam_value :
             max_beam_value = beams['XX'][xx,yy]
             max_beam_value_xx = xx
             max_beam_value_yy = yy


    #Get Haslam and interpolate onto grid
    my_map=get_Haslam(frequency,radio_image=radio_image)      
    mask=np.isnan(za_grid)
    za_grid[np.isnan(za_grid)]=90.0 #Replace nans as they break the interpolation     
    sky_grid=map_sky_astropy(my_map['skymap'],my_map['RA'],my_map['dec'],gps,az_grid,za_grid,epoch='B1950') # B1950 due to HASLAM being in this one !
    sky_grid[mask]=np.nan #Remask beyond the horizon
    ra_map = my_map['RA']
    dec_map = my_map['dec']
    
    if max_beam_value_xx >= 0 and max_beam_value_yy >= 0 :
       print("Maximum beam value %.8f is at (xx,yy) = (%d,%d) -> (az,za) = (%.4f,%.4f) , T_sky = %.2f [K]" % (max_beam_value,max_beam_value_xx,max_beam_value_yy,phi[max_beam_value_xx,max_beam_value_yy],theta[max_beam_value_xx,max_beam_value_yy],sky_grid[max_beam_value_xx,max_beam_value_yy]))
       # ERROR : IndexError: too many indices       print "Maximum beam value %.8f is at (xx,yy) = (%d,%d) -> (az,za) = (%.4f,%.4f) , T_sky = %.2f [K] at (RA,DEC) = (%.4f,%.4f) [deg]" % (max_beam_value,max_beam_value_xx,max_beam_value_yy,phi[max_beam_value_xx,max_beam_value_yy],theta[max_beam_value_xx,max_beam_value_yy],sky_grid[max_beam_value_xx,max_beam_value_yy],ra_map[max_beam_value_xx,max_beam_value_yy],dec_map[max_beam_value_xx,max_beam_value_yy])
       
    

    beamsky_sum_XX=0
    beam_sum_XX=0
    Tant_XX=0    
    beam_dOMEGA_XX=0
    beam_dOMEGA_sum_XX=0
    beamsky_sum_YY=0
    beam_sum_YY=0
    Tant_YY=0    
    beam_dOMEGA_YY=0    
    beam_dOMEGA_sum_YY=0                            

    out_file=open(out_filename,"a+")
    for pol in pols:
        # Get gridded sky
        print('frequency',frequency)
        beam=beams[pol]
        beamsky=beam*sky_grid
        beam_dOMEGA=beam*dOMEGA
        
        beam_dOMEGA_sum=np.nansum(beam_dOMEGA)
        beamsky_sum=np.nansum(beamsky)
        beam_sum=np.nansum(beam)
        Tant=np.nansum(beamsky)/np.nansum(beam)

        print('sum(beam)',np.nansum(beam))
        print('sum(beamsky)',np.nansum(beamsky))                                        
        print('Tant=sum(beamsky)/sum(beam)=', Tant)

        if pol == 'XX' :
           beamsky_sum_XX = beamsky_sum
           beam_sum_XX    = beam_sum   
           Tant_XX        = Tant
           beam_dOMEGA_sum_XX = beam_dOMEGA_sum
                                               
        if pol == 'YY' :
           beamsky_sum_YY = beamsky_sum
           beam_sum_YY    = beam_sum   
           Tant_YY        = Tant
           beam_dOMEGA_sum_YY = beam_dOMEGA_sum                                                                                                    

        # add receiver temperature :
        if use_trcv :
           print("T_ant = %.2f K ( T_ant_XX = %.2f K , T_ant_YY = %.2f K), T_rcv = %.2f" % (Tant,Tant_XX,Tant_YY,T_rcv))
           Tant_XX = Tant_XX + T_rcv*(1.00 + rcv_noise_coupling ) # plus noise with some coupling        
           Tant_YY = Tant_YY + T_rcv*(1.00 + rcv_noise_coupling ) # plus noise with some coupling
           Tant = Tant + T_rcv*(1.00 + rcv_noise_coupling ) # plus noise with some coupling

           print("After + T_rcv and coupling : T_ant = %.2f K ( T_ant_XX = %.2f K , T_ant_YY = %.2f K )" % (Tant,Tant_XX,Tant_YY))
           

        uxtime=gps+315964800
        out_file.write('%d %.8f %.8f\n' % (uxtime,Tant,lst))
        
        
        filename='%s_%.2fMHz_%s_%s' % (gps,frequency/1.0e6,pol,model)
        fstring="%.2f"%(frequency/1.0e6)
            
        
        if plottype=='all':
            plottypes=['beam','sky','beamsky','beamsky_scaled'] 
        else:
            if plottype is None or plottype.lower() =='none' :
                plottypes=[]
            else :
                plottypes=[plottype]            
            
        for pt in plottypes:
            if pt =='beamsky':
                textlabel='Beam x sky %s (LST %.2f hr), %s MHz, %s-pol, Tant=%.1f K' % (gps, get_LST(gps),fstring, pol, Tant)
                plot_beamsky(beamsky,frequency, textlabel, filename, extension,
                             figsize=figsize,directory=directory,lst=lst)
            elif pt =='beamsky_scaled':        
                textlabel='Beam x sky (scaled) %s (LST %.2f hr), %s MHz, %s-pol, Tant=%.1f K (max T=%.1f K)' % \
                            (gps, get_LST(gps),fstring, pol, Tant,np.nanmax(beamsky))                
                plot_beamsky(beamsky,frequency, textlabel, filename+'_scaled', extension, 
                         figsize=figsize, vmax=np.nanmax(beamsky)*0.4,directory=directory,lst=lst)    
    
            elif pt =='beam':
                textlabel='Beam for %s, %s MHz, %s-pol' % (gps,fstring, pol)
                plot_beamsky(beam,frequency, textlabel, filename+'_beam', extension,
                             figsize=figsize, cbar_label='' ,directory=directory, lst=lst)
            elif pt =='sky':
                textlabel='Sky for %s (LST %.2f hr), %s MHz, %s-pol' % (gps, get_LST(gps),fstring, pol)
                plot_beamsky(sky_grid, frequency, textlabel, filename+'_sky', extension,
                             figsize=figsize,directory=directory, lst=lst)
    
    out_file.close()
    
    if False :
       # This is integration test to verify if this way of integrating beam is better than using projection and if it agrees with showspec 
       # see memo : /home/msok/Desktop/AAVS2/logbook/20200915_array_factor_beam_vs_feko_beam.odt
       test_sum_domega = 0
       test_sum = 0
       test_sum2 = 0
       # integration test :
       theta_deg = 0
       step = 1 # 0.1 on bighorns
       step_rad = step*(math.pi/180.00)
       while theta_deg <= 90.00 :
          print("DEBUG_INTEGRATION : theta_deg = %.4f [deg] , test_sum = %.10f , test_domega = %.10f ( test_sum2 = %.10f )" % (theta_deg,test_sum,test_sum_domega,test_sum2))
          theta_rad = (math.pi/180.00)*theta_deg
          phi_deg = 0.00       
          while phi_deg <= 360.00 :
             phi_rad = (math.pi/180.00)*phi_deg          
             dOmega = math.sin(theta_rad)*step_rad*step_rad
          
             beams_xx,beams_yy=eda_beam.get_eda_beam( numpy.array([[theta_rad]]), numpy.array([[phi_rad]]), resolution=resolution, zenithnorm=True, power=True, freq=frequency, lst=lst, 
                                                      pointing_za_deg=pointing_za_deg, pointing_az_deg=pointing_az_deg, gain_sigma_dB=gain_sigma_dB, 
                                                      gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, dipole_type=dipole_type, xpos=xpos, ypos=ypos, zpos=zpos,
                                                      use_beam_fits=use_beam_fits, station_name=station_name, projection=projection )

             test_sum += beams_yy[0,0]*dOmega
             test_sum2 += beams_yy[0,0]
             test_sum_domega += dOmega
          
             phi_deg += step
          
          theta_deg += step
       
       print("DEBUG_INTEGRATION : test_sum = %.10f , test_domega = %.10f ( test_sum2 = %.10f )" % (test_sum,test_sum_domega,test_sum2))
    
    return (beamsky_sum_XX,beam_sum_XX,Tant_XX,beam_dOMEGA_sum_XX,beamsky_sum_YY,beam_sum_YY,Tant_YY,beam_dOMEGA_sum_YY,beams)


#FIXME: most the arguments in make_primarybeammap are not needed
def get_beam_power(  gps, delays, frequency, model, pointing_az_deg=0, pointing_za_deg=0, zenithnorm=True, gain_sigma_dB=0.0, gain_sigma_ph_160mhz=0.00, dipole_type='short', xpos=None, ypos=None, zpos=None, 
                     use_beam_fits=False, station_name="EDA", projection="zea" ):

   lst=get_LST(gps)   
   # first go from altitude to zenith angle
   theta_rad=pointing_za_deg*math.pi/180   
   phi_rad=pointing_az_deg*math.pi/180     
               
   theta=np.array([[theta_rad]])
   phi=np.array([[phi_rad]])                       

   beams={}   
   beams['XX'],beams['YY']=eda_beam.get_eda_beam( theta, phi, zenithnorm=True, power=True, freq=frequency, lst=lst, pointing_za_deg=pointing_za_deg, pointing_az_deg=pointing_az_deg, 
                                                  gain_sigma_dB=gain_sigma_dB, gain_sigma_ph_160mhz=gain_sigma_ph_160mhz, dipole_type=dipole_type, xpos=xpos, ypos=ypos, zpos=zpos,
                                                  use_beam_fits=use_beam_fits, station_name=station_name, projection=projection )
   
   return beams
       
    

def plot_beamsky(beamsky, frequency, textlabel, filename, extension,
                 figsize=8, vmax=None, cbar_label='beam x Tsky (K)', 
                 directory=None, dec=-26.7033, lst=0, 
                 projection='ZEA' ):       
    # do the plotting
    # this sets up the figure with the right aspect ratio
    
    fig=pylab.figure(figsize=(figsize,0.6*figsize),dpi=300)
    pylab.axis('on')
    ax1=fig.add_subplot(1,1,1,polar=False)
    
    pylab.axis('off')
    # Add polar grid on top (but transparent background)
#TODO: change grid labels to ZA.
    ax2=fig.add_subplot(1,1,1,polar=True, frameon=False)
    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.patch.set_alpha(0.0)    
    ax2.tick_params(color='0.5', labelcolor='0.5')
    for spine in ax2.spines.values():
        spine.set_edgecolor('0.5')        
    ax2.grid(which='major', color='0.5') 
    
    #Beamsky example:
    if vmax is not None:
        im=ax1.imshow(beamsky,interpolation='none', vmax=vmax)
    else:
        im=ax1.imshow(beamsky,interpolation='none')
    #Add colorbar on own axis
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(im, cax=cbar_ax,label='Tsky (K)') #Require more recent numpy (e.g. 1.9.2 works)
    cbar=fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(cbar_label)

    ax1.set_title(textlabel+'\n\n')

    full_filename=filename
    if directory is not None:
       full_filename=directory + '/' + filename
    try:
        fig.savefig(full_filename+'.'+extension) #transparent=True if we  want transparent png 
    except RuntimeError as err:
        logger.error('Error saving figure: %s\n' % err)
        return None

    # save fits files:
    full_filename=filename + '.fits'
    if directory is not None:
       full_filename=directory + '/' + filename + '.fits'
    print("Filename2 = %s" % filename)        
    try:
        hdu = pyfits.PrimaryHDU()
#        beamsky(np.isnan(hdu.data))=0.00
        hdu.data = beamsky

        # nan -> 0
        hdu.data[np.isnan(hdu.data)]=0.00
        
        # add keywords:
#        pixscale=180.0/(beamsky.shape[0]/2) # for all-sky
# based on eda_beam.py :

# 2020-10-15 - this will only work ok for SIN (slant orhtographic correction ) not ZEA, but ZEA is required for equal area pixels in the integration 
#        dec_rad = math.fabs( -26.98138889*(math.pi/180.00) ) # WARNING : this one uses zenith MRO DEC at EPOCH =  B1950 as HASLAM is in this epoch and otherwise there is shift by a few pixels !
        dec_rad = math.fabs( dec*(math.pi/180.00) ) # WARNING : Haslam epoch B1950 vs. J2000 - still to be sorted !!!
        pixscale = 180.0/ math.pi
        mu = 0.00 #  1e20 - not ok # 0.00 ~ ok, 1.00 - not ok
        # see Eq. 16 in Greisen and Carballeta /home/msok/Desktop/SKA/Tim/doc/AllSky_imaging/references
#        pixscale_final = (180.0/ (math.pi*beamsky.shape[0])) * ( (mu+1) * math.cos(dec_rad) ) / ( mu + 	math.sin(dec_rad) )
        if projection == 'ZEA' :
           pixscale_final = (180.0/ (math.pi*((beamsky.shape[0])/math.sqrt(2.00)))) * ( (mu+1) * math.cos(dec_rad) ) / ( mu +         math.sin(dec_rad) )
        else :
           pixscale_final = (180.0/ (math.pi*beamsky.shape[0])) * ( (mu+1) * math.cos(dec_rad) ) / ( mu + math.sin(dec_rad) )  

        CTYPE1_string = 'RA---%s' % (projection)
        CTYPE2_string = 'DEC---%s' % (projection)
                
        dec_rad = math.fabs( dec*(math.pi/180.00) )
        hdu.header['CRPIX1'] = beamsky.shape[0]/2 # + 1
        hdu.header['CDELT1'] = pixscale_final # see Eq. 16 in Greisen and Carballeta /home/msok/Desktop/SKA/Tim/doc/AllSky_imaging/references
        hdu.header['CRVAL1'] = lst*15.00 
        hdu.header['CTYPE1'] = CTYPE1_string # or SIN for SIN projection see 348 call in (az_grid, za_grid,n_total,dOMEGA) = fits_beam.makeAZZA(resolution,projection='ZEA')
        hdu.header['CRPIX2'] = beamsky.shape[0]/2 # + 1
        hdu.header['CDELT2'] = pixscale_final # see Eq. 16 in Greisen and Carballeta /home/msok/Desktop/SKA/Tim/doc/AllSky_imaging/references
        hdu.header['CRVAL2']  = dec
        hdu.header['CTYPE2']  = CTYPE2_string # or SIN for SIN projection see 348 call in (az_grid, za_grid,n_total,dOMEGA) = fits_beam.makeAZZA(resolution,projection='ZEA')          
        hdu.header['BEAM_AZ'] = 0 
        hdu.header['BEAM_ZA'] = 0 
        hdu.header['FREQ']    = frequency
        hdu.header['EPOCH']   = 2000.00 # 1950.00 # 2020-10-15 fix : This is because HASLAM is in B1950 and otherwise sources have position shifted by a few pixels !
                       
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(full_filename,clobber=True)        
        print("Saved output image to file %s" % full_filename)
    except RuntimeError as err:
        logger.error('Error saving figure: %s\n' % err)
        return None
   

    pylab.close()
    
######################################################################
# Running as executable
if __name__=='__main__':
    main()
