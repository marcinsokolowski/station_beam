from django import forms
from django.utils.safestring import mark_safe

class ParametersForm(forms.Form):
    lst_hours = forms.FloatField( label='LST (0 - 24) [hours]' , min_value=0.0 , max_value=24.0, initial=0.00 ) # ,max_digits=5,decimal_places=1)
    azimuth_deg   = forms.FloatField( label='Azimuth (0 - 360) [deg]'   , min_value=0.0 , max_value=360.0, initial=0.00 ) # ,max_digits=5,decimal_places=1)
    elevation_deg = forms.FloatField( label='Elevation (0 - 90) [deg]' , min_value=0.0 , max_value=90.0, initial = 90.00 ) # ,max_digits=5,decimal_places=1)

    STATION_CHOICES = [
       ('EDA2', 'EDA2'),
       ('AAVS2', 'AAVS2')
    ] 
    station_name = forms.CharField( label='Station' , widget=forms.Select(choices=STATION_CHOICES) , initial='eda2' )   
    
    CHOICES = [('1', 'Show plot'), ('2', 'Save zip file (data and plots)')]    
    mode = forms.ChoiceField( widget=forms.RadioSelect , choices=CHOICES , label='Output format' , initial='1' )
