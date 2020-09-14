from __future__ import print_function
import beam_tools
# from pylab import *
import numpy
import math
import sys

import astropy.io.fits as pyfits
from pylab import loadtxt 

from scipy.interpolate import griddata
import scipy.interpolate
import matplotlib.pyplot as plt


read_feko_debug_level = 0

max_theta = -1000
max_phi   = -1000

def angular_distance_degree( phi1_deg, theta1_deg, phi2_deg, theta2_deg ) :
   phi1_rad = phi1_deg*(math.pi/180.00)
   phi2_rad = phi2_deg*(math.pi/180.00)
   theta1_rad = theta1_deg*(math.pi/180.00)
   theta2_rad = theta2_deg*(math.pi/180.00)
   
   dist_rad_cos = math.sin(theta1_rad)*math.sin(theta2_rad) + math.cos(theta1_rad)*math.cos(theta2_rad)*math.cos( phi1_rad - phi2_rad )
   dist_rad = math.acos( dist_rad_cos )
   
   return (dist_rad*(180.00/math.pi))



def read_feko_file( feko_file="AAVS1_full_station_160MHz_WithoutTurnedOFF-ANTS.txt" ) :
    global max_phi
    global max_theta

    theta_phi = loadtxt( feko_file , usecols = (0,1), comments=["#","*"] )
    values    = loadtxt( feko_file , usecols = (2,) , comments=["#","*"])
    data      = loadtxt( feko_file , comments=["#","*"] )
    
    theta_deg = theta_phi[:,0]
    phi_deg   = theta_phi[:,1]

    theta_unique = numpy.unique(theta_deg)
    phi_unique   = numpy.unique(phi_deg)
    
    max_theta = max( theta_unique )
    max_phi   = max( phi_unique )
    
    ( theta_mesh, phi_mesh ) = numpy.meshgrid(theta_unique, phi_unique)
    beam_values_2d = theta_mesh.copy()*0.00
    
    step_theta = theta_unique[1] - theta_unique[0]
    step_phi   = phi_unique[1]   - phi_unique[0]
    
    print("read_feko_file(%s) : step_theta = %.4f, step_phi = %.4f, beam_values_2d.shape = %s" % (feko_file,step_theta,step_phi,beam_values_2d.shape))

    max_value = -1000
    max_t_deg = -1000
    max_p_deg = -1000    
    for i in range(0,len(values)) :
        t_deg = theta_deg[i]
        p_deg = phi_deg[i]
        beam   = values[i]
        
        theta_idx = int( round( t_deg / step_theta ) )
        phi_idx   = int( round( p_deg / step_phi ) )
        
        beam_values_2d[phi_idx,theta_idx] = beam 
        
        if beam > max_value :
           max_value = beam
           max_t_deg = t_deg
           max_p_deg = p_deg
        
        print("\t\t%d : (%.4f,%.4f) = %.4f -> theta_idx = %d, phi_idx = %d" % (i,t_deg,p_deg,beam,theta_idx,phi_idx))
        
    print("Max value = %.2f at (theta,phi) = (%.4f,%.4f) [deg] , Angular range : max_phi = %.4f [deg], max_theta = %.4f [deg]" % (max_value,max_t_deg,max_p_deg,max_phi,max_theta))    
        
    return (theta_phi,values,data, theta_unique, phi_unique, step_theta, step_phi, theta_mesh, phi_mesh, beam_values_2d)

def find_beam_value( phi_deg, th_deg, step_theta, step_phi, beam_values_2d, interpolate=False, interpolation_radius=2, theta_map=None, phi_map=None, interpolate_to_nearest=False ) :
   if step_theta == 0 or step_phi == 0 :
      print("ERROR : step_theta = %.4f , step_phi = %.4f" % (step_theta,step_phi))
      return -10000

   if numpy.isnan(th_deg) :
      print("ERROR : th_deg = %.4f" % (th_deg))
      return -10000

   if numpy.isnan(phi_deg) :
      print("ERROR : phi_deg = %.4f" % (phi_deg))
      return -10000

   th_index  = int( round( th_deg / step_theta ) )
   phi_index = int( round( phi_deg / step_phi ) )
   
   if phi_index >= beam_values_2d.shape[0] :
      if phi_deg > max_phi and math.fabs( phi_deg - 360.00 ) < 0.0000001 :
         print("DEBUG : phi_deg=%.2f > max_phi = %.2f [deg] -> phi_deg := phi_index := 0" % (phi_deg,max_phi))
         phi_deg = 0.00
         phi_index = 0
         
      if phi_index == 360 and beam_values_2d.shape[0] == 360 :
         print("DEBUG : phi_deg=%.2f , phi_index=%d , max_phi = %.2f [deg] -> phi_deg := phi_index := 0" % (phi_deg,phi_index,max_phi))
         phi_index = 0
         phi_deg = 0.00
   
#   if theta_map is not None and phi_map is not None :
#       print "Map shapes = %s and %s" % (theta_map.shape,phi_map.shape)

   beam = 0.00   
   if phi_index >= 0 and phi_index < beam_values_2d.shape[0] and th_index >= 0 and th_index < beam_values_2d.shape[1] :                
       beam = beam_values_2d[phi_index,th_index]
       
       if interpolate :
          # interpolation of beam value using pixels in interpolation_radius :
          if theta_map is not None or phi_map is not None :
              beam_sum = 0.00
              weigth_sum = 0.00
              min_dist_deg = 100000.00
              beam_at_nearest = 0.00
              
#              x_list=[]
#              y_list=[]
#              z_list=[]
              
              continue_sum = True
              for yy in range(th_index-interpolation_radius,th_index+interpolation_radius+1) :
                 for xx in range(phi_index-interpolation_radius,phi_index+interpolation_radius+1) :
                     if xx>=0 and xx<beam_values_2d.shape[0] and yy>=0 and yy<beam_values_2d.shape[1] :
#                         if yy>=30 :
#                             print "DEBUG : (x,y) = (%d,%d)" % (xx,yy)
                         phi_xx = (phi_map[xx,yy])*(180.00/math.pi)
                         th_xx  = (theta_map[xx,yy])*(180.00/math.pi)
                         
#                         x_list.append(xx)
#                         y_list.append(yy)
#                         z_list.append( beam_values_2d[xx,yy] )
    
                     
                         dist_deg = angular_distance_degree( phi_deg, th_deg, phi_xx, th_xx )                
                         
                         # if beam value found in a distance < 1 arcmin, use it :
                         if (dist_deg) < (1.00/60.00) :
                             beam_sum = beam_values_2d[xx,yy]
                             weigth_sum = 1.00
                             continue_sum = False
                             
                         if continue_sum : 
                             beam_sum   = beam_sum + beam_values_2d[xx,yy] / dist_deg
                             weigth_sum = weigth_sum + 1.00 / dist_deg
                         
                         if dist_deg < min_dist_deg :
                            min_dist_deg = dist_deg
                            beam_at_nearest = beam_values_2d[xx,yy]
                            
              if interpolate_to_nearest :
                  beam = beam_at_nearest
              else : 
                  # works but won't be easy to convert from (th_deg,phi_deg) -> (th_deg_x_value, phi_deg_x_value ) for interpolation
                  # points = numpy.array([x_list,y_list])
                  # values = numpy.array(z_list)
                  # interpolator = scipy.interpolate.CloughTocher2DInterpolator( points , values )
                  # beam = interpolator( th_deg_x_value, phi_deg_x_value )
                                    
                  # weigthed by distance            
                  beam = beam_sum / weigth_sum
              
              
          
          else :
              print("ERROR : cannot use interpolation of beam as theta and phi maps are not provided (=None)")       
           
       
       if read_feko_debug_level > 0 :
           print("DEBUG :  find_beam_value( %.4f, %.4f, %.4f, %.4f , beam_values_2d ) : th_index = %d, phi_index = %d, beam = %.8f" % (phi_deg, th_deg, step_theta, step_phi, th_index, phi_index, beam))
       
       
       if beam > 1 :
          print("WARNING : wrong beam value = %.4f at (phi_deg,theta_deg) = (%.4f,%.4f) [deg] , step_phi = %.2f, step_theta = %.2f , phi_index = %d, theta_index = %d" % (beam,phi_deg,th_deg,step_phi,step_theta,phi_index,th_index))
   else :
       print("ERROR in find_beam_value( %.4f, %.4f, %.4f, %.4f , beam_values_2d ) -> phi_index = %d or th_index = %d outside of range" % (phi_deg, th_deg, step_theta, step_phi, phi_index, th_index))
       beam = -10000       
   
   return beam

def find_beam_value_slowest( az_value, th_value, azim_rad_array, theta_rad_array, beam_values ) :
    l = len( azim_rad_array )
    
    ret = -100000
    min_dist = 1e20
    for i in range(0,l) :
       az = azim_rad_array[i]
       th = theta_rad_array[i]
    
       dist_rad_cos = math.sin(th)*math.sin(th_value) + math.cos(th)*math.cos(th_value)*math.cos( az - az_value )
       dist_rad = math.acos( dist_rad_cos )

       if dist_rad < min_dist :
          min_dist = dist_rad
          ret = beam_values[i]
          
    return ret      
           

# theta,phi in radians ! 
def get_feko_beam( theta, phi, feko_file="AAVS1_full_station_160MHz_WithoutTurnedOFF-ANTS.txt", interpolate=False, interpolation_radius=1, interpolate_to_nearest=False ) :    
    beam = numpy.zeros( theta.shape )
    
    # read beam form file :
    (theta_phi,values,data, theta_unique, phi_unique, step_theta, step_phi, theta_mesh, phi_mesh, beam_values_2d) = read_feko_file( feko_file=feko_file )    
    
    # fill all beam pixels :           
    for y in range(0,phi.shape[0]) :
       print("PROGRESS : y = %d / %d" % (y,phi.shape[0]))
       
       for x in range(0,phi.shape[1]) :              
          az_deg = phi[x,y]*(180.00/math.pi)
          th_deg = theta[x,y]*(180.00/math.pi)
          
          # convert az_deg to phi_feko_deg
          phi_feko_deg = 90.00 - az_deg
          if phi_feko_deg < 0 :
             phi_feko_deg = phi_feko_deg + 360.00
          if phi_feko_deg > 360.00 :
             phi_feko_deg = phi_feko_deg - 360.00
          
          if numpy.isnan(th_deg) :
              if read_feko_debug_level > 0 :
                  print("ERROR at (x,y) = (%d,%d), th_deg = NaN, phi_deg = %.2f [deg]" % (x,y,phi_feko_deg))
              beam[x,y] = 0.00
          else :          
              beam_value = find_beam_value( phi_feko_deg, th_deg, step_theta, step_phi, beam_values_2d, theta_map=theta_mesh, phi_map=phi_mesh, interpolate=interpolate, interpolation_radius=interpolation_radius, interpolate_to_nearest=interpolate_to_nearest )
              beam[x,y] = beam_value
              
#          print "DEBUG : beam[%d,%d] = %.4f at az_deg = %.2f -> phi_deg = %.2f , th_deg = %.2f" % (x,y,beam[x,y],az_deg,phi_feko_deg,th_deg)
                            
    
    return (beam)

def get_feko_beam_slow( theta, phi, feko_file="AAVS1_full_station_160MHz_WithoutTurnedOFF-ANTS.txt" ) :    
    beam = numpy.zeros( theta.shape )
    
    # read beam form file :
    (theta_phi,values,data) = read_feko_file( feko_file=feko_file )    
    
    theta_deg = theta_phi[:,0]
    phi_deg   = theta_phi[:,1]
    
    print("Size of theta_deg = %d, phi_deg = %d" % (len(theta_deg),len(phi_deg)))
    
    azim_deg = 90.00 - phi_deg

    # convert phi_feko -> azim     
    for i in range(0,azim_deg.shape[0]) :
        if azim_deg[i] < 0 :
           azim_deg[i] = azim_deg[i] + 360.00
        if azim_deg[i] > 360.00 :
           azim_deg[i] = azim_deg[i] - 360.00

    azim_rad  = azim_deg*(math.pi/180.00)
    theta_rad = theta_deg*(math.pi/180.00)

    # fill all beam pixels :           
    for y in range(0,phi.shape[0]) :
       print("y = %d / %d" % (y,phi.shape[0]))
       for x in range(0,phi.shape[1]) :              
          az_rad = phi[x,y]
          th_rad = theta[x,y]
          
          beam_value = find_beam_value( az_rad, th_rad, azim_rad, theta_rad, values, theta_map=theta, phi_map=phi )
          beam[x,y] = beam_value
                            
    
    return (beam,azim_rad,theta_rad)


    
def main() :
   feko_file = "AAVS2_256_elem_station_160MHz_Xpol.txt"
   if len(sys.argv) > 1:
       feko_file = sys.argv[1]

   (theta_phi,values,data, theta_unique, phi_unique, step_theta, step_phi, theta_mesh, phi_mesh, beam_values_2d) = read_feko_file( feko_file=feko_file )
   size_x = beam_values_2d.shape[0]
   size_y = beam_values_2d.shape[1]
   
   hdu = pyfits.PrimaryHDU()
   hdu.data = beam_values_2d
   hdulist = pyfits.HDUList([hdu])
   hdulist.writeto('new.fits',clobber=True)
   
   # plt.subplot(221)
   # plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
   # plt.plot(phi_theta, values, 'k.', ms=1)
   # plt.title('Original')   
   # plt.show()
    
if __name__ == "__main__":
    main()
    