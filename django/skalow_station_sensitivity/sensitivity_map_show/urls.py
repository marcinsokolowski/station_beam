from django.urls import path
from .views import sensitivity_map_show

urlpatterns = [ 
    path('sensitivity_map_show',sensitivity_map_show),
]       

