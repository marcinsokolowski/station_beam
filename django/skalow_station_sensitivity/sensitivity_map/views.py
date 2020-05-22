from django.shortcuts import render
from .models import Post
from .forms import ParametersForm

# Create your views here.

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

    return render(request, 'parameters_map.html', {'form': form})
