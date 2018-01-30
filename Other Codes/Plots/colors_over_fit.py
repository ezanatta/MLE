import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

fits_file = 'fits-null.fits'
hdu = fits.open(fits_file)[0]
wcs = WCS(hdu.header)

RA = np.loadtxt('2768mag.dat', usecols=(2,))
DEC = np.loadtxt('2768mag.dat', usecols=(3,))

with open('2768GC.dat', 'r') as f:
    RAstr = [line.split()[0] for line in f]
with open('2768GC.dat', 'r') as f:    
    DECstr = [line.split()[1] for line in f]

Rc = np.loadtxt('2768mag.dat', usecols=(4,))
Rcerr = np.loadtxt('2768mag.dat', usecols=(5,))
z = np.loadtxt('2768mag.dat', usecols=(6,))
zerr = np.loadtxt('2768mag.dat', usecols=(7,))

colour = Rc - z
c_sep = 0.57   #colour separation
n = len(RA)

RAblue = []
DECblue= []
RAred = []
DECred = []

for i in range(0,n):
    if colour[i] < c_sep:
        RAblue.append(RA[i])
        DECblue.append(DEC[i])
    else:
        RAred.append(RA[i])
        DECred.append(DEC[i])

coords = zip(RA, DEC)
coords_blue = zip(RAblue, DECblue)
coords_red = zip(RAred, DECred)

ra = []
dec = []
ra_b = []
dec_b = []
ra_r = []
dec_r = []


pix = wcs.wcs_world2pix(coords, 1)
pix_b = wcs.wcs_world2pix(coords_blue, 1)
pix_r = wcs.wcs_world2pix(coords_red, 1)

for i in range(0,len(pix)):
    ra.append(pix[i][0])
    dec.append(pix[i][1])
for i in range(0,len(pix_b)):
    ra_b.append(pix_b[i][0])
    dec_b.append(pix_b[i][1])
for i in range(0, len(pix_r)):
    ra_r.append(pix_r[i][0])
    dec_r.append(pix_r[i][1])




plt.plot(ra_r, dec_r, 'ro')
plt.plot(ra_b, dec_b, 'bo')
plt.imshow(hdu.data,cmap='gray_r', label='All GCs')
plt.show()
