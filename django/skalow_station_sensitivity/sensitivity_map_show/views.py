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
def sensitivity_map_show(request):
   zipfile = False
   post = Post.objects.all()

   params = None   
   if request.method == 'GET':
      print("GET !") 
      params = request.GET

   elif request.method == 'POST':
      print("POST !!!")
      params = request.POST

   lst_hours = float( params['lst_hours'] )
   frequency_mhz = float( params['frequency_mhz'] )
   mode          = int( params['mode'] )
   return_zip_file = False
   if mode >= 2 :
       return_zip_file = True
   print("DEBUG : station = %s" % (params['station_name']))
   (station) = params['station_name'] # "EDA2"
   print("DEBUG : station = %s" % (station))
   db_path = ( "%s/" % (config.sensitivity_db_path) )
   
   print("Parameters = %s -> LST = %.4f [hours], frequency = %.2f [MHz] , station = %s, db_path = %s , mode = %d ( return_zip_file = %s)" % (params,lst_hours,frequency_mhz,station,db_path,mode,return_zip_file))
   
   out_fitsname_base = "sensitivity_map_lst%06.2fh_freq%06.2fMHz" % (lst_hours, frequency_mhz)
   ( azim_x , za_x , aot_x , sefd_x , azim_y , za_y , aot_y , sefd_y, out_txt_filename_X, out_txt_filename_Y ) = sensitivity_db.get_sensitivity_map( 
                                                                                                                frequency_mhz, lst_hours, output_file_base=out_fitsname_base, db_path=db_path, 
                                                                                                                station=station, out_fitsname_base=out_fitsname_base, output_dir=config.save_output_path )


#   return render(request,"sensitivity_vs_freq_show/index.html"  )

   # 
   save_output_path = config.save_output_path
   mkdir_p( save_output_path )

   # create plots :   
   print("DEBUG : calling sensitivity_db.plot_sensitivity_map (saving to %s)" % (save_output_path))
   output_file_base = out_fitsname_base
#   (png_image_path,buf) = sensitivity_db.plot_sensitivity( freq_x, aot_x, freq_y, aot_y, output_file_base=output_file_base )
#   def plot_sensitivity_map( azim_deg, za_deg, aot, out_fitsname_base="sensitivity" , save_text_file=False, do_plot=False, freq_mhz=0.00, lst_h=0.00, pol="Unknown" ) :
   ( lut_x , sensitivity_x, sens_fits_file_x , azim_fits_file_x , za_fits_file_x , out_txt_filename_x , pngfile_x, buf_x )  = sensitivity_db.plot_sensitivity_map( azim_x, za_x, aot_x, out_fitsname_base=output_file_base, save_text_file=False, do_plot=True, freq_mhz=frequency_mhz, lst_h=lst_hours, pol="X", out_dir=config.save_output_path )
   ( lut_y , sensitivity_y, sens_fits_file_y , azim_fits_file_y , za_fits_file_y , out_txt_filename_y , pngfile_y, buf_y )  = sensitivity_db.plot_sensitivity_map( azim_y, za_y, aot_y, out_fitsname_base=output_file_base, save_text_file=False, do_plot=True, freq_mhz=frequency_mhz, lst_h=lst_hours, pol="Y", out_dir=config.save_output_path )

#   return render(request,"sensitivity_vs_freq_show/index.html"  )
   
   
#   out_file_name = save_output_path + "/" + output_file_base
#   (text_file) = sensitivity_db.save_sens_vs_freq_file( freq_x, aot_x, sefd_x, freq_y, aot_y, sefd_y, out_file_name )   


#   return render(request,"sensitivity_vs_freq_show/index.html" )
   
   buf_x.seek(0)
   string_x = base64.b64encode(buf_x.read())
   uri_x = 'data:image/png;base64,' + urllib.parse.quote(string_x)

   buf_y.seek(0)
   string_y = base64.b64encode(buf_y.read())
   uri_y = 'data:image/png;base64,' + urllib.parse.quote(string_y)

   response = None
   zip_file_name = None
   if zipfile :
      zip_file_name = save_output_path + "/" + output_file_base + ".zip"
      in_memory = StringIO()
      zip = ZipFile( zip_file_name , "w" )
      print("DEBUG : created zip archive %s" % (zip_file_name))
      print("DEBUG : adding file %s to the zip archive as %s" % (text_file , os.path.basename(text_file) ))
      
      zip.write( out_txt_filename_X , os.path.basename(out_txt_filename_X) )
      zip.write( out_txt_filename_Y , os.path.basename(out_txt_filename_Y) )
      zip.write( pngfile_x , os.path.basename(pngfile_x)  )
      zip.write( pngfile_y , os.path.basename(pngfile_y)  )

      # fix for Linux zip files read in Windows
      for file in zip.filelist:
         file.create_system = 0

      zip.close()
      
      if return_zip_file : # working - if needed :
         response = FileResponse(open( zip_file_name , 'rb'))

         return response

   lst_hours_str = "%.2f" % (lst_hours)
   freq_str      = "%.2f" % (frequency_mhz)
   args = { 'image_x':uri_x , 'image_y':uri_y , 'zipfile':zip_file_name, 'lst':lst_hours_str, 'frequency':freq_str }
   print("DEBUG : mode = %d" % (mode))
         
   return render(request,"sensitivity_map_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   render(request,"sensitivity_vs_lst_show/index.html" , args ) # , context_instance=RequestContext(request) )
#   return response
