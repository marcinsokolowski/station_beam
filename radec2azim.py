import sys
import math

from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)

DEG2RAD = (math.pi/180.00)
RAD2DEG = (180.00/math.pi)

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
   
   return (azim_deg,alt_deg)

# https://docs.astropy.org/en/stable/coordinates/index.html
def radec2gal( ra_deg, dec_deg ) :
   c_icrs = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
   
   return (c_icrs.galactic.l.value,c_icrs.galactic.b.value)

def gal2radec( glon_deg, glat_deg ) :
   c_icrs = SkyCoord(l=glon_deg*u.degree, b=glat_deg*u.degree, frame='galactic')
   
   return (c_icrs.icrs.ra.value,c_icrs.icrs.dec.value)
   
   

def main() :
   ra_deg = 0.00
   if len(sys.argv) > 1:
       ra_deg = float( sys.argv[1] )

   dec_deg = 0.00
   if len(sys.argv) > 2:
       dec_deg = float( sys.argv[2] )

   lst_hours = 0.00
   if len(sys.argv) > 3:
       lst_hours = float( sys.argv[3] )
       
   geo_lat = MWA_POS.lat.value
   if len(sys.argv) > 4:
      geo_lat = float( sys.argv[4] )
      
   radec2azim( ra_deg , dec_deg , lst_hours , geo_lat )   
      
if __name__ == "__main__":
    main()
