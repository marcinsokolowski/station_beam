import numpy
import math

def is_digit( znak ) :
   if znak=='0' or znak=='1' or znak=='2' or znak=='3' or znak=='4' or znak=='5' or znak=='6' or znak=='7' or znak=='8' or znak=='9' :
      return True

   return False


def round_and_get_error( value , rel_error=0.2, debug=False ):
   error = rel_error*value 
   error_string = ("%.20f" % error)
   x = error
   
   n=1   
   error_string = ("%.20f" % x)
   len_buff = len(error_string)
   first_digit = -1
   after_dot=-1
   
   for i in range(0,len_buff): 
      if is_digit( error_string[i] ) :
         if after_dot >= 0 :
            after_dot += 1

         if first_digit < 0 and error_string[i]!='0' :
            first_digit = error_string[i];
            break;
      else :
         if error_string[i] == '.' :
            after_dot = 0
      
   rounded = x;

#   after_dot++;
# COMMENT IFS for a)
   if after_dot < 0 :
      after_dot = 0

   if first_digit != '1' :
      # rounded = ceil(x*TMath::Power(10,after_dot))/TMath::Power(10,after_dot);
      rounded = math.ceil(x*math.pow(10,after_dot))/(math.pow(10,after_dot))

   if debug :
      print("First digit in %s is %c, after dot = %d -> rounded = %.20f\n" % (error_string,first_digit,after_dot,rounded))
   if first_digit == '1' :
      n = 2

#   double scale = TMath::Power(10.0, ceil(log10(fabs(x))) + n);
#     return round(x * scale) / scale;
   if after_dot > 0 :
       error_string = ( "%.*g" % (n, rounded))
   else :
       error_string = ("%.0f" % rounded)

   if debug :
      print("Buff = %s\n" % error_string)

   dot_idx=error_string.find(".")
   if dot_idx >= 0 :
      digits=error_string[dot_idx+1:]
      if debug :
         print("DEBUG : error_string := %s , digits = %s" % (error_string,digits))
      len_digits=len(digits)
   else :
      len_digits=0

   # len_digits=len(error_string)
   ret = ('%.*f' % (len_digits,round(value,len_digits)))

   return (ret,error_string)
      
   


def round_and_get_error_OLD( value , rel_error=0.2, debug=False ):
   error = rel_error*value        
   error_string = ("%.20f" % error)
   
   if error > 1 :
      error_string = str(round(error))
   else :  
      error_len=len(error_string)
      last_digit=-1
      for i in range(0,error_len):
         if error_string[i] != '0' and error_string[i] != '.' :
#            if error_string[i] == '1' :
            last_digit=i
            break
   
      if debug : 
         print("DEBUG : error_string = %s -> last_digit = %d" % (error_string,last_digit))
         
      if last_digit >= 0 :
         if error_string[last_digit] == '1' :               
            last_digit += 1
            while last_digit < error_len and error_string[last_digit]=='0' :
               last_digit += 1
         
         error_string = error_string[0:last_digit+1]

   error_len=len(error_string)
   dot_idx=error_string.find(".")
   if dot_idx >= 0 :
      digits=error_string[dot_idx+1:]
      if debug :
         print("DEBUG : error_string := %s , length = %d , digits = %s" % (error_string,error_len,digits))
      len_digits=len(digits)
   else :
      len_digits=0

   # ret = ('%.*g' % (error_len,value))
   ret = ('%.*f' % (len_digits,round(value,len_digits)))
   
#   ret = ('%.20f' % value)
#   ret = ('%f' % round(value,n_digits))
         
   return (ret,error_string)
   
        