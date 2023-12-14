from django.urls import path, include, re_path

from . import views

urlpatterns = [
    path('', views.main, name='main'),
    re_path(r'autocomplete/', views.autocomplete_view, name='autocomplete_view'),
    re_path(r'results/(.+)/(.+)', views.results, name='results'),
    re_path(r'results/(.+)', views.results, name='results'),    
]
