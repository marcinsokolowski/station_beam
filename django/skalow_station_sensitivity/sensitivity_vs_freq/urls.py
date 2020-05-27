from django.urls import path
from .views import get_parameters

urlpatterns = [ 
#    path('sensitivity_vs_lst',blog_list),
    path('sensitivity_vs_freq',get_parameters),
]       
