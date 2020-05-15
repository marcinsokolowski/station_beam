from django.urls import path
from .views import sensitivity_vs_lst_show

urlpatterns = [ 
#    path('sensitivity_vs_lst',blog_list),
    path('sensitivity_vs_lst_show',sensitivity_vs_lst_show),
]	
