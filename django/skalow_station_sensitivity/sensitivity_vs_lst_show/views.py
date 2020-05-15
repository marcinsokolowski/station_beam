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

   print("Parameters = %s" % params)      

   return render(request,"sensitivity_vs_lst_show/index.html")
