import sys
import os
import errno
from string import Template

# from django.shortcuts import render
from .models import Post
from django.shortcuts import *
from django.template import RequestContext
from django.http import FileResponse

sys.path.append("../")       # for config
# sys.path.append("../../../") # for sensitivity_db
import config
sys.path.append( config.sensitivity_db_path )
print("DEBUG : added path %s" % (config.sensitivity_db_path))
sys.path.append( config.station_beam_path )
print("DEBUG : added path %s" % (config.station_beam_path))

# do not require DISPLAY :
import matplotlib
if 'matplotlib.backends' not in sys.modules:
    matplotlib.use('agg')

# from io import BytesIO
import sensitivity_db
sensitivity_db.init_web_interface()

# import cStringIO
from zipfile import ZipFile
from io import BytesIO
from io import StringIO
import io
import urllib, base64

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise


# Create your views here.
def sensitivity_radec_astro_show(request):
   zipfile = False
   post = Post.objects.all()

   params = None   
   if request.method == 'GET':
      print("GET !") 
      params = request.GET

   elif request.method == 'POST':
      print("POST !!!")
      params = request.POST

   frequency_mhz = float( params['frequency_mhz'] )
   ra_deg   = float( params['ra_deg'] )
   dec_deg = float( params['dec_deg'] )
   ha_start_h = float( params['ha_start_h'] )
   ha_end_h = float( params['ha_end_h'] )
   bw_mhz = float( params['bw_mhz'] )
   n_stations = int( params['n_stations'] )
   inttime = float( params['inttime'] )
   mode          = int( params['mode'] )
   plot_type = params['plot_type']
   return_zip_file = False
   if mode >= 2 :
       zipfile = True
       return_zip_file = True
   print("DEBUG : station = %s" % (params['station_name']))
   (station) = params['station_name'] # "EDA2"
   print("DEBUG : station = %s" % (station))
   db_path = ( "%s/" % (config.sensitivity_db_path) )
   
   print("Parameters = %s -> %.4f MHz, (ra,dec) = (%.4f,%.4f) [deg] , HA range = %.4f - %.4f [h], BW = %.4f [MHz], Inttime = %.2f [s], N_stations = %d, station = %s, db_path = %s , mode = %d ( return_zip_file = %s)" %  (params,frequency_mhz,ra_deg,dec_deg,ha_start_h,ha_end_h,bw_mhz,inttime,n_stations,station,db_path,mode,return_zip_file))
   
   lst_start = ha_start_h + ra_deg/15.00
   lst_end = ha_end_h + ra_deg/15.00
   
   if ha_start_h <= -12 and ha_end_h >= +12 :
      lst_start = 0.00
      lst_end   = 24.00

   if lst_start < 0 :
      lst_start = 24.00 + lst_start 
   
   if lst_end < 0 :
      lst_end = 24.00 + lst_end
   
   if lst_start <= lst_end :
      print("DEBUG : OK lst_start = %.8f [h] < lst_end = %.8f [h]" % (lst_start,lst_end))
   else :
      # lst_end < lst_start -> adding 24 hours to lst_end, but it has to be dealt later inside this function below :
      lst_end += 24 
      
   print("DEBUG : sensitivity_radec_astro_show LST range = %.8f - %.8f [hours]" % (lst_start,lst_end))
   
   # TODO : add a single number inside this functions !!!   
   (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y, lst_i,aot_i,sefd_i, noise_x, noise_y, noise_i, noise_x_total, noise_y_total, noise_i_total,total_integration_time ) = sensitivity_db.get_sensitivity_radec_lstrange( ra_deg, dec_deg, frequency_mhz, lst_start=lst_start, lst_end=lst_end, time_step=inttime, station=station , db_path=db_path, bandwidth_hz=bw_mhz*1e6, n_stations=n_stations )

   
   # 
   save_output_path = config.save_output_path
   mkdir_p( save_output_path )

   noise_info = None
   if noise_x_total < 1.00 and noise_y_total < 1.00 and noise_i_total < 1.00 :   
      noise_info = "Expected %.4f [hours] image noise for X/Y/Stokes I is  %.6f [mJy] , %.6f [mJy] , %.6f [mJy]" % (total_integration_time,noise_x_total*1000.00,noise_y_total*1000.00,noise_i_total*1000.00)
   else :
      noise_info = "Expected %.4f [hours] image noise for X/Y/Stokes I is  %.3f [Jy] , %.3f [Jy] , %.3f [Jy]" % (total_integration_time,noise_x_total,noise_y_total,noise_i_total)

   # create plots :   
   if plot_type == 'LST' :
      print("DEBUG : calling sensitivity_db.plot_sensitivity_vs_lst (saving to %s)" % (save_output_path))
      output_file_base = "%s_sensitivity_lst0-24h_ra%.2fdeg_dec%.2fdeg_%.2fMHz" % (station,ra_deg,dec_deg,frequency_mhz)
      (png_image_path,buf) = sensitivity_db.plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, lst_start=lst_start, lst_end=lst_end, azim_deg=-1000, za_deg=-1000, freq_mhz=frequency_mhz, output_file_base=output_file_base, do_show=False, save_output_path=save_output_path, lst_i=lst_i, aot_i=aot_i, ra_deg=ra_deg, dec_deg=dec_deg, station_name=station, noise_info=noise_info )
   else : 
      print("DEBUG : calling sensitivity_db.plot_sensitivity_vs_ha (saving to %s)" % (save_output_path))
      output_file_base = "%s_sensitivity_vs_ha_ra%.2fdeg_dec%.2fdeg_%.2fMHz" % (station,ra_deg,dec_deg,frequency_mhz)
      (png_image_path,buf) = sensitivity_db.plot_sensitivity_vs_ha( lst_x, aot_x, lst_y, aot_y, ha_start_h=ha_start_h, ha_end_h=ha_end_h, azim_deg=-1000, za_deg=-1000, freq_mhz=frequency_mhz, output_file_base=output_file_base, do_show=False, save_output_path=save_output_path, lst_i=lst_i, aot_i=aot_i, ra_deg=ra_deg, dec_deg=dec_deg, station_name=station, noise_info=noise_info )

   
   # def save_sens_vs_lst_file( lst_x, aot_x, sefd_x, lst_y, aot_y, sefd_y out_file_base ) :
   out_file_name = save_output_path + "/" + output_file_base
   (text_file) = sensitivity_db.save_sens_vs_lst_file( lst_x, aot_x, sefd_x, lst_y, aot_y, sefd_y, out_file_name )   

#   print("DEBUG : plotting image %s" % (png_image_path))
   
#   template = loader.get_template('sensitivity_vs_lst_show/index.html')

#   return render_to_response('sensitivity_vs_lst_show/index.html', { 'image': png_image_path }, context_instance=RequestContext(request) )
   
#   {'form': form,'obs': observation, 
#                                                           'setting': setting, 'error': error,
#                                                           'image': image_path, 'obs_link': obs_link}, 
#                                                           context_instance=RequestContext(request))

   buf.seek(0)
   string = base64.b64encode(buf.read())
   uri = 'data:image/png;base64,' + urllib.parse.quote(string)

   response = None
   zip_file_name = None
   if zipfile :
      zip_file_name = save_output_path + "/" + output_file_base + ".zip"
      in_memory = StringIO()
#      zip = ZipFile(in_memory, "w")
      zip = ZipFile( zip_file_name , "w" )
      print("DEBUG : created zip archive %s" % (zip_file_name))
      print("DEBUG : adding file %s to the zip archive as %s" % (text_file , os.path.basename(text_file) ))
      
#      zipinfo = zip.ZipInfo('helloworld.txt', date_time=time.localtime(time.time()))
#      zipinfo.create_system = 1
#      zip.writestr( "test.txt" , StringIO('helloworld').getvalue())
#      zip.writestr(  StringIO('helloworld').getvalue() , "test.txt" )
#      zip.writestr( "test.txt" , StringIO('helloworld') )
      
#      io.open(text_file) 
#      zip.writestr( text_file, os.path.basename(text_file) )

      zip.write( text_file , os.path.basename(text_file) )
      zip.write( png_image_path , os.path.basename(png_image_path)  )
#      zip.writestr( "test.jpg" , uri  )

      # flagged_tiles_file = "%d_flagged_tiles.txt" % cal_id
      # zip.writestr( flagged_tiles_file , flagged_tiles )

      # fix for Linux zip files read in Windows
      for file in zip.filelist:
         file.create_system = 0

      zip.close()
      
      if return_zip_file : # working - if needed :
         response = FileResponse(open( zip_file_name , 'rb'))
#        response = HttpResponse( content_type="application/zip" )
#        response["Content-Disposition"] = "attachment; filename=%s" % ( output_file_base + ".zip" )
#        in_memory.seek(0)
  #      response.write(in_memory.read())

         return response


#   print("DEBUG string = %s -> %s" % (string,uri))
#    args = {'form':form, 'text':text, 'image':uri}
   # if noise_x_total < 1.00 and noise_y_total < 1.00 and noise_i_total < 1.00 :
   
   noise_x_str = ""
   noise_y_str = ""
   noise_i_str = ""
   if noise_x_total < 1.00 and noise_y_total < 1.00 and noise_i_total < 1.00 :
      noise_x_str = "%.6f [mJy]" % (noise_x_total*1000.00)
      noise_y_str = "%.6f [mJy]" % (noise_y_total*1000.00)
      noise_i_str = "%.6f [mJy]" % (noise_i_total*1000.00)
   else :
      noise_x_str = "%.6f [Jy]" % (noise_x_total)
      noise_y_str = "%.6f [Jy]" % (noise_y_total)
      noise_i_str = "%.6f [Jy]" % (noise_i_total)
   
   args = { 'image':uri , 'zipfile':zip_file_name, 'noise_x_str':noise_x_str, 'noise_y_str':noise_y_str, 'noise_i_str':noise_i_str, 'warning':"" }
   print("DEBUG : mode = %d" % (mode))
   
   if bw_mhz > 20 :
      args['warning'] = ('WARNING : bandwidth = %.2f MHz but sensitivity at the center frequency used' % (bw_mhz))   

   return render(request,"sensitivity_radec_astro_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   render(request,"sensitivity_vs_lst_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   return response