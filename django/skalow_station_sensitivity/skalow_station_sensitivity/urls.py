"""skalow_station_sensitivity URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
# from sensitvity_vs_lst import views

from . import views

urlpatterns = [
#    path('', views.index, name='index'),
    path('admin/', admin.site.urls),
    path('', include('hello_world.urls')),
    path('', include('sensitivity_vs_lst.urls')),
    path('', include('sensitivity_vs_lst_show.urls')),
    path('', include('sensitivity_vs_freq.urls')),
    path('', include('sensitivity_vs_freq_show.urls')),
    path('', include('sensitivity_radec_vs_freq.urls')),
    path('', include('sensitivity_radec_vs_freq_show.urls')),
    path('', include('sensitivity_map.urls')),
    path('', include('sensitivity_map_show.urls')),
]
