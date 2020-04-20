from django.urls import path
from .views import blog_list

urlpatterns = [ 
   path('sensitivity_vs_lst',blog_list)
]	
