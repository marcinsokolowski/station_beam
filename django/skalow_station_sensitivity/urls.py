from django.urls import path

from . import views

urlpatterns = [
    path('posts/', sensitivity_vs_lst.urls),
    path('posts/', sensitivity_vs_lst_show.urls)
]
