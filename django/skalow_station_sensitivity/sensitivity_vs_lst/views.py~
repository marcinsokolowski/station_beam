from django.shortcuts import render
from .models import Post

from .forms import ParametersForm

def get_name(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = NameForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            return HttpResponseRedirect('/thanks/')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = NameForm()

    return render(request, 'name.html', {'form': form})

def get_parameters(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = ParametersForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            return HttpResponseRedirect('/thanks/')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = ParametersForm()

    return render(request, 'parameters.html', {'form': form})

# Create your views here.
def blog_list(request):
   post = Post.objects.all()
   context = {
      'blog_list':post
   }
   
   return render(request,"sensitivity_vs_lst/index.html",context)


def hello_world(request):
    return render(request, 'sensitivity_vs_lst/hello_world.html', {})

   