from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import sys
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import galaxy_inputs


#this code plots all tracers, colored, on a fits file provided. 
#                INPUTS and INFO:
#----> Up to three catalogs: stars, PNe and GCs, and a fits file ( one big enough to fit every tracer inside).
#
#----> Most of time the catalogs I used had different formats and units. 
#      The PN.S catalogs are in this fashion: hh:mm:ss. So I need to load them as strings,
#      use SkyCoord to load these strings as RA and DEC and correct units do decimal degrees.
#----->GC catalogs from Vicenzo Pota are already in decimal degrees, so I load them directly with numpy.
#----->Star catalogs will differ (I guess) between each other, but one of this two methods will work for any case.
#-----> Of course you can use Topcat or any other way to correct the coordinates also. 
#       The important thing is to have everything in decimal degrees so astropy.fits can understand it.

def read_catalog_gc(gccat):
    RAg = np.loadtxt(gccat, usecols=(1,))
    DECg = np.loadtxt(gccat, usecols=(2,))
    
    return RAg, DECg
    
def read_catalog_pne(pnecat):
    RAp = np.loadtxt(pnecat, usecols=(1,))
    DECp = np.loadtxt(pnecat, usecols=(2,))
    
    return RAp, DECp
    
def read_catalog_stars(starcat):
    with open(starcat, 'r') as f:
        RAs = [line.split()[0] for line in f]               
    with open(starcat, 'r') as f:    
        DECs = [line.split()[1] for line in f]
    
    for i in range(0, len(RAs)):
	c = SkyCoord(RAs[i], DECs[i])
	RAs[i], DECs[i] = c.ra.degree, c.dec.degree
    
    return RAs, DECs    


gal = raw_input('Enter the galaxy number: ')
galinput = galaxy_inputs(gal)
gc = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'GC.dat'  
pne = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'PNE.dat'  
stars = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'STARS.dat'

op = raw_input('Do we have GC? y/n: ')
op2 = raw_input('Do we have PNe? y/n: ')
op3 = raw_input('Do we have Stars? y/n: ')

if op == 'y':
    RA, DEC = read_catalog_gc(gc)
    coords = zip(RA, DEC)
    ra = []
    dec = []
if op2 == 'y':
    RApne, DECpne = read_catalog_pne(pne)
    coordsPNE = zip(RApne, DECpne)
    ra_pne = []
    dec_pne = []
if op3 == 'y':
    RAstar, DECstar = read_catalog_stars(stars)
    coords_star = zip(RAstar, DECstar)
    ra_star = []
    dec_star = []

fits_file = '/home/emilio/MLE/Galaxies/'+gal+'/'+'fits-null.fits'
hdu = fits.open(fits_file)[0]
wcs = WCS(hdu.header)

if op == 'y':
    pix = wcs.wcs_world2pix(coords, 1)    
if op3 =='y':                                                    #here we transform the RA and DEC coords in pixel coords 
    pix_star = wcs.wcs_world2pix(coords_star, 1)                #to overplot on the fits file in the right places. 
if op2 == 'y':
    pix_pne = wcs.wcs_world2pix(coordsPNE, 1)                   #This is done using the WCS loaded before.
                                                            #so it is important that the fits you used had the WCS.

#if you really need to use a fits file without WCS (like when you plot the f-map), check how it is done in gen_cat.py
if op == 'y':
    for i in range(0,len(pix)):
        ra.append(pix[i][0])
        dec.append(pix[i][1])
if op3 =='y': 
    for i in range(0,len(pix_star)):
        ra_star.append(pix_star[i][0])
        dec_star.append(pix_star[i][1])
if op2 == 'y':
    for i in range(0, len(pix_pne)):
        ra_pne.append(pix_pne[i][0])
        dec_pne.append(pix_pne[i][1])    

if op == 'y' and op2 == 'y' and op3 == 'y':
    allgc = plt.subplot(1,1,1, projection=wcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r') #you can change the contrast in the cmap options.
    plt.plot(ra, dec, marker='o',markerfacecolor='None', linestyle='none', markeredgewidth=1, markeredgecolor='magenta', label='GCs')
    plt.plot(ra_star, dec_star, marker='s', color='orange', markeredgewidth=1, markeredgecolor='orange',markerfacecolor='None', linestyle='none', label='Stars')
    plt.plot(ra_pne, dec_pne, marker='o', markeredgewidth=1,markeredgecolor='green',markerfacecolor='None',linestyle='none', label='PNe')
    plt.legend(loc='lower right')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()

if op == 'y' and op2 == 'y' and op3 != 'y':
    allgc = plt.subplot(1,1,1, projection=wcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r') 
    plt.plot(ra, dec, marker='o',markerfacecolor='None', linestyle='none', markeredgewidth=3, markeredgecolor='magenta', label='GCs')
    plt.plot(ra_pne, dec_pne, marker='o', markeredgewidth=2,markeredgecolor='green',markerfacecolor='None',linestyle='none', label='PNe')
    plt.legend(loc='lower right')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()
    
if op2 == 'y' and op != 'y' and op3 != 'y':
    allgc = plt.subplot(1,1,1, projection=wcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r') 
    plt.plot(ra_pne, dec_pne, marker='o', markeredgewidth=2,markeredgecolor='green',markerfacecolor='None',linestyle='none', label='PNe')
    plt.legend(loc='lower right')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()
    
if op == 'y' and op2 != 'y' and op3 != 'y':
    allgc = plt.subplot(1,1,1, projection=wcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r') #you can change the contrast in the cmap options.
    plt.plot(ra, dec, marker='o',markerfacecolor='None', linestyle='none', markeredgewidth=1, markeredgecolor='magenta', label='GCs')
    plt.legend(loc='lower right')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()
