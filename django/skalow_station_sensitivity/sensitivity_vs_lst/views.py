from django.shortcuts import render
from .models import Post

# Create your views here.
def blog_list(request):
   post = Post.objects.all()
   context = {
      'blog_list':post
   }
   
   return render(request,"sensitivity_vs_lst/index.html",context)


def hello_world(request):
    return render(request, 'sensitivity_vs_lst/hello_world.html', {})

   