from django.shortcuts import render
from .models import Post


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
   
   print("Parameters = %s -> %.4f MHz, (az,el) = (%.4f,%.4f) [deg]" % (params,frequency_mhz,azimuth_deg,elevation_deg))

   return render(request,"sensitivity_vs_lst_show/index.html")
