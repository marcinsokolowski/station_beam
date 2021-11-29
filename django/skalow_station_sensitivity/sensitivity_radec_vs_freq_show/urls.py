from django.urls import path
from .views import sensitivity_radec_vs_freq_show

urlpatterns = [ 
#    path('sensitivity_vs_lst',blog_list),
    path('sensitivity_radec_vs_freq_show',sensitivity_radec_vs_freq_show),
]       

