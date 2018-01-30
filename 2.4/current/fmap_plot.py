import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

DRA = np.loadtxt('2MASS.GC.cat', usecols=(0,))
DDEC = np.loadtxt('2MASS.GC.cat', usecols=(1,))
#f = np.loadtxt('N2768_f.GC.dat', usecols=(1,))


fits_file = 'f.fits'
hdu = fits.open(fits_file)[0]
wcs = WCS(hdu.header)

coords = zip(DRA, DDEC)

ra = []
dec = []

pix = wcs.wcs_world2pix(coords, 1)

for i in range(0,len(pix)):
    ra.append(pix[i][0])
    dec.append(pix[i][1])


plt.imshow(hdu.data, origin='lower',cmap='gray', label='All GCs')
#plt.scatter(DRA, DDEC, c=f, cmap='Blues')
plt.plot(DRA, DDEC, 'bo')
#cb = plt.colorbar()
#cb.set_label('Probability Of Being in the Spheroid')
plt.title('PNE')
plt.show()





