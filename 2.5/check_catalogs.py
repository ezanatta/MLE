import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


Vari = np.loadtxt('NGC2768.GC.red.ML.dat', usecols=(2, ))
#Vemilio = np.loadtxt('red_catalogue.dat', usecols=(3,))

RAari = np.loadtxt('NGC2768.GC.red.ML.dat', usecols=(0,))
DECari = np.loadtxt('NGC2768.GC.red.ML.dat', usecols=(1,))

RAemilio = np.loadtxt('red_catalogue.dat', usecols=(0, ))
DECemilio = np.loadtxt('red_catalogue.dat', usecols=(1, ))

RAemilio = RAemilio*u.degree
DECemilio = DECemilio*u.degree

RAemilio = RAemilio.to(u.arcsec)
DECemilio = DECemilio.to(u.arcsec)
RAemilio = RAemilio.value
DECemilio = DECemilio.value


RAari = np.around(RAari)
RAemilio = np.around(RAemilio)

for i in range(0, len(RAemilio)):
	c = RAemilio[i]
	g = RAari[i]
	if g != c:
		print i
		print 'my gc ra', c
		print 'arianna gc ra', g

 
