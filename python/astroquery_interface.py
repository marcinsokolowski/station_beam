
has_astroquery = False
try :
   from astroquery.simbad import Simbad
   from astropy.coordinates import SkyCoord
   from astropy.coordinates import ICRS, Galactic, FK4, FK5
   import astropy.units as u
   
   has_astroquery = True
except :
   print("ERROR : something wrong with import or call to astroquery")
   has_astroquery = False


def get_radec_for_a_source( source_name ) :
   if has_astroquery :
      ra_deg = None
      dec_deg = None
      found = False
   
      try :
         source_info = Simbad.query_object( source_name )
         found = True
         # general :
         # source_info = Simbad.query_object( source_name )
      except :
         print("WARNING : could not find object %s in SIMBAD catalogue" % (source_name))
         found = False
      
      if found and source_info is not None :
         coord_str=source_info['RA'][0] + " " + source_info['DEC'][0]
         coordinates = SkyCoord(coord_str, frame=FK5, unit=(u.hourangle, u.deg), obstime="J2000")

         ra_deg=coordinates.ra.value
         dec_deg=coordinates.dec.value

         print("DEBUG : Right ascension = %.4f [deg] , declination = %.4f [deg]" % (ra_deg,dec_deg))
      else : 
         print("WARNING : source %s not found in the catalogue" % source_name)
      
      return (ra_deg,dec_deg)
   else :
      print("WARNING : astroquery is not available (install with pip install astroquery)")

   return (None,None)      
   