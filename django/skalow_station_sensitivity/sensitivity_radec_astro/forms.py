from django import forms
from django.utils.safestring import mark_safe

class NameForm(forms.Form):
    your_name = forms.CharField(label='Your name', max_length=100)
    
    
class ParametersForm(forms.Form):
    frequency_mhz = forms.FloatField( label='Center frequency [MHz]' , min_value=40.0 , max_value=350.0, initial=160.00 ) # ,max_digits=5,decimal_places=1)
    ra_deg   = forms.FloatField( label='RA [deg]'   , min_value=0.0 , max_value=360.0, initial=0.00 ) # ,max_digits=5,decimal_places=1)
    dec_deg = forms.FloatField( label='DEC [deg]' , min_value=-90.0 , max_value=90.0, initial = 0.00 ) # ,max_digits=5,decimal_places=1)
    ha_start_h = forms.FloatField( label='Hour Angle Start [hours]'   , min_value=-12 , max_value=0.0, initial=-2.00 ) 
    ha_end_h = forms.FloatField( label='Hour Angle End [hours]'   , min_value=0, max_value=12.0, initial=2.00 ) 
    bw_mhz = forms.FloatField( label='Observing Bandwidth [MHz]'   , min_value=0.00, max_value=300, initial=30.00 )
    n_stations = forms.IntegerField( label='Number of stations'   , min_value=512 , max_value=512, initial=512 )
    inttime = forms.FloatField( label='Integrations [seconds]'   , min_value=0.5, max_value=3600, initial=120 ) # duration of individual integrations 

    STATION_CHOICES = [
       ('EDA2', 'EDA2'),
       ('AAVS2', 'AAVS2')
    ] 
    station_name = forms.CharField( label='Station' , widget=forms.Select(choices=STATION_CHOICES) , initial='eda2' )   
    
    PLOT_CHOICES = [
       ('Hour_Angle', 'Hour_Angle'),
       ('LST', 'LST')
    ] 
    plot_type = forms.CharField( label='Plot A/T [m^2/K] vs.' , widget=forms.Select(choices=PLOT_CHOICES) , initial='Hour_Angle' )
    
    CHOICES = [('1', 'Show image'), ('2', 'Save zip file (data and image)')]    
    mode = forms.ChoiceField( widget=forms.RadioSelect , choices=CHOICES , label='Output format' , initial='1' )
