import sys
import os
import errno

# from django.shortcuts import render
from .models import Post
from django.shortcuts import *
from django.template import RequestContext
from django.http import FileResponse

sys.path.append("../")
import config
sys.path.append( config.station_beam_path )

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
def sensitivity_vs_lst_show(request):
   zipfile = True
   post = Post.objects.all()

   params = None   
   if request.method == 'GET':
      print("GET !") 
      params = request.GET

   elif request.method == 'POST':
      print("POST !!!")
      params = request.POST

   frequency_mhz = float( params['frequency_mhz'] )
   azimuth_deg   = float( params['azimuth_deg'] )
   elevation_deg = float( params['elevation_deg'] )
   mode          = int( params['mode'] )
   return_zip_file = False
   if mode >= 2 :
       return_zip_file = True
   za_deg = (90.00 - elevation_deg)
   station = "EDA2"
   db_path = ( "%s/" % (config.sensitivity_db_path) )
   
   print("Parameters = %s -> %.4f MHz, (az,el) = (%.4f,%.4f) [deg] , station = %s, db_path = %s , mode = %d ( return_zip_file = %s)" % (params,frequency_mhz,azimuth_deg,elevation_deg,station,db_path,mode,return_zip_file))
   
   (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y) = sensitivity_db.get_sensitivity_lstrange( azimuth_deg, za_deg, frequency_mhz, lst_start=0, lst_end=24, time_step=300, station="EDA2" , db_path=db_path )

   
   # 
   save_output_path = config.save_output_path
   mkdir_p( save_output_path )

   # create plots :   
   print("DEBUG : calling sensitivity_db.plot_sensitivity_vs_lst (saving to %s)" % (save_output_path))
   output_file_base = "%s_sensitivity_lst0-24h_az%.2fdeg_za_%.2fdeg_%.2fMHz" % (station,azimuth_deg,za_deg,frequency_mhz)
   (png_image_path,buf) = sensitivity_db.plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, lst_start=0, lst_end=20, azim_deg=azimuth_deg, za_deg=za_deg, freq_mhz=frequency_mhz, output_file_base=output_file_base, do_show=False, save_output_path=save_output_path )
   
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
   args = { 'image':uri , 'zipfile':zip_file_name }
   print("DEBUG : mode = %d" % (mode))
         
   return render(request,"sensitivity_vs_lst_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   render(request,"sensitivity_vs_lst_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   return response
