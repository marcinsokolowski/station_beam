import sys

from django.shortcuts import render
from .models import Post

sys.path.append("../")
import config
sys.path.append( config.station_beam_path )
import sensitivity_db


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

   return render(request,"sensitivity_vs_lst_show/index.html")
