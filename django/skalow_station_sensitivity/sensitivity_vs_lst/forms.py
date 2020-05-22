from django import forms
from django.utils.safestring import mark_safe

class NameForm(forms.Form):
    your_name = forms.CharField(label='Your name', max_length=100)
    
    
class ParametersForm(forms.Form):
    frequency_mhz = forms.FloatField(label='Frequency [MHz]',min_value=0.0, max_value=350.0) # ,max_digits=5,decimal_places=1)
    azimuth_deg   = forms.FloatField(label='Azimuth [deg]',min_value=0.0, max_value=360.0) # ,max_digits=5,decimal_places=1)
    elevation_deg = forms.FloatField(label='Elevation [deg]',min_value=0.0, max_value=90.0) # ,max_digits=5,decimal_places=1)

    STATION_CHOICES = [
       ('eda2', 'EDA2'),
       ('aavs2', 'AAVS2')
    ] 
    station_name = forms.CharField(label='Station', widget=forms.Select(choices=STATION_CHOICES))   
    
    CHOICES = [('1', 'Show image'), ('2', 'Save zip file (data and image)')]    
    mode = forms.ChoiceField(widget=forms.RadioSelect, choices=CHOICES, label='mode', initial='1')
