import sys
from astropy import time

MJD = float(sys.argv[1])
t = time.Time(MJD, format='mjd')
print int(t.unix)

