from django.urls import path
from .views import blog_list
from .views import get_parameters

urlpatterns = [ 
#    path('sensitivity_vs_lst',blog_list),
    path('sensitivity_radec_astro',get_parameters),
]	
