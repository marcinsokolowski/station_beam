from django import forms

class NameForm(forms.Form):
    your_name = forms.CharField(label='Your name', max_length=100)
    
    
class ParametersForm(forms.Form):
    frequency_mhz = forms.FloatField(label='Frequency [MHz]',min_value=0.0, max_value=350.0) # ,max_digits=5,decimal_places=1)
    
    