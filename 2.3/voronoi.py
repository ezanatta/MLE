#import matplotlib
#matplotlib.use('Agg')               #if remotely connect decomment this
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt

from cap_display_bins_generators import display_bins_generators

filename = "Cappellari2011a_Atlas3D_Paper1_UnBinned_Cubes_v1.0/MS_NGC5866_r7_C2D.fits"
hdu = pyfits.open(filename)
spectrum = hdu[0].data
table = hdu[2].data

x = table["A"] # Coordinates of the original spaxels in arcsec (nort is up)
y = table["D"]
flux = np.mean(spectrum, 1) # surface brightness

filename = "atlas3d_stellar_kinematics/cappellari2011a/PXF_bin_MS_NGC5866_r7_idl.fits.gz"
hdu = pyfits.open(filename)
table = hdu[1].data

xgen = table["XS"] # Voronoi generators
ygen = table["YS"]
velbin = table["VPXF"] # Mean stellar velocity

display_bins_generators(xgen, ygen, velbin, x, y)
plt.tricontour(x, y, -2.5*np.log10(flux/np.max(flux)), levels=np.arange(20)) # 1 mag contours
plt.show()