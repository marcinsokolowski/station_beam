import sys

# from django.shortcuts import render
from .models import Post
from django.shortcuts import *
from django.template import RequestContext

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


# Create your views here.
def sensitivity_vs_lst_show(request):
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
   za_deg = (90.00 - elevation_deg)
   station = "EDA2"
   db_path = ( "%s/" % (config.sensitivity_db_path) )
   
   print("Parameters = %s -> %.4f MHz, (az,el) = (%.4f,%.4f) [deg] , station = %s, db_path = %s" % (params,frequency_mhz,azimuth_deg,elevation_deg,station,db_path))
   
   (lst_x,aot_x,sefd_x, lst_y,aot_y,sefd_y) = sensitivity_db.get_sensitivity_lstrange( azimuth_deg, za_deg, frequency_mhz, lst_start=0, lst_end=24, time_step=300, station="EDA2" , db_path=db_path )

   
   # 
   save_output_path = config.save_output_path

   # create plots :   
   print("DEBUG : calling sensitivity_db.plot_sensitivity_vs_lst (saving to %s)" % (save_output_path))
   output_file_base = "%s_sensitivity_lst0-24h_az%.2fdeg_za_%.2fdeg_%.2fMHz" % (station,azimuth_deg,za_deg,frequency_mhz)
   (png_image_path) = sensitivity_db.plot_sensitivity_vs_lst( lst_x, aot_x, lst_y, aot_y, lst_start=0, lst_end=20, azim_deg=azimuth_deg, za_deg=za_deg, freq_mhz=frequency_mhz, output_file_base=output_file_base, do_show=False, save_output_path=save_output_path )

   print("DEBUG : plotting image %s" % (png_image_path))
   
#   template = loader.get_template('sensitivity_vs_lst_show/index.html')

#   return render_to_response('sensitivity_vs_lst_show/index.html', { 'image': png_image_path }, context_instance=RequestContext(request) )
   
#   {'form': form,'obs': observation, 
#                                                           'setting': setting, 'error': error,
#                                                           'image': image_path, 'obs_link': obs_link}, 
#                                                           context_instance=RequestContext(request))


   return render(request,"sensitivity_vs_lst_show/index.html" , { 'image': png_image_path } ) # , context_instance=RequestContext(request) )
