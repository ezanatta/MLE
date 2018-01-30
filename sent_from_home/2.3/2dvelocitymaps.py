from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import galaxy_inputs
from termcolor import colored
import rpy2.robjects as robjects

#to do: eliminate the y axis on second and third plots

robjects.r('''
    plotVor <- function(){
    
    library(deldir)
    library(datautils)
    
    gal <- paste(readLines('~/MLE/2.2/current_galaxy.dat'), collapse=" ")
    path_to_cat = paste0('~/MLE/Galaxies/',gal, '/')
    cat_name = paste0(path_to_cat, '/N',gal,'GC.dat')
    
    cat = read.table(cat_name)
    
    RA = cat$V2
    DEC = cat$V3
    V = cat$V4
    g = cat$V5
    i = cat$V6
    color = g-i
    
    df = data.frame(RA, DEC)
    voronoi = deldir(df$RA, df$DEC)
    rbPal <- colorRampPalette(c('blue','red'))
    df$Col <- rbPal(10)[as.numeric(cut(V,breaks = 10))]
    
    X11()
    
    plot.deldir(voronoi, fill=df$Col, wlines='tess', wpoints='none', xlab='RA', ylab='DEC', main=paste0('NGC', gal))
    dev.copy(png,filename=paste0(path_to_cat,'plots/voronoi', gal,'.png'))
    }
''')

def read_catalog_gc(galcat):
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    
    V = np.loadtxt(galcat, usecols=(3,))        
    
    col1 = np.loadtxt(galcat, usecols=(4,))
    col2 = np.loadtxt(galcat, usecols=(5,))
    colour = col1 - col2
    
    return RA, DEC, V, colour  
    
def read_catalog_pne(galcat):
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    
    V = np.loadtxt(galcat, usecols=(3,))        
    
    colour = -1000
    
    return RA, DEC, V, colour 
    
plotVor = robjects.r['plotVor']

gal = raw_input('Enter the galaxy number: ')
galinput = galaxy_inputs(gal)
choice = raw_input('GC or PNE? ')
if choice == 'GC':
    galcat = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'GC.dat'   

    with open('current_galaxy.dat', 'w+') as current:     
        print >>current, gal     

    RAg, DECg, V, col = read_catalog_gc(galcat)
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]    
    
    c_sep = float(inp[4])   #colour separation
    n = len(RAg)
else:
    galcat = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'PNE.dat'   

    with open('current_galaxy.dat', 'w+') as current:     
        print >>current, gal     
        
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]

    RAg, DECg, V, col = read_catalog_pne(galcat)
    
    c_sep = col   #colour separation
    n = len(RAg)

if c_sep == -1000:       
    coords = zip(RAg, DECg)
    
    fits_file = '/home/emilio/MLE/Galaxies/'+gal+'/'+'f.fits'
    hdu = fits.open(fits_file)[0]
    wcs = WCS(hdu.header)
    
    ra = []
    dec = []
    
    pix = wcs.wcs_world2pix(coords, 1)
    
    for i in range(0,len(pix)):
        ra.append(pix[i][0])
        dec.append(pix[i][1])
    
    Vmax = max(V)-800
    Vmin = min(V)
    
    fig = plt.subplot(1,1,1, projection=wcs)
    plt.scatter(ra, dec, c=V, vmin=Vmin, vmax=Vmax)
    cb = plt.colorbar(fraction=0.04, orientation='horizontal')
    cb.set_label('Velocity $(km/s)$')
    plt.imshow(hdu.data,cmap='gray_r', label='All GCs')
    plt.title('All GCs')
    
    fig.set_ylabel('RA')
    fig.set_xlabel('DEC')
    
    plt.show()
    
else:
    RAblue = []
    DECblue= []
    RAred = []
    DECred = []
    V_b = []
    V_r = []
    
    for i in range(0,n):
        if col[i] < c_sep:
            RAblue.append(RAg[i])
            DECblue.append(DECg[i])
            V_b.append(V[i])
        else:
            RAred.append(RAg[i])
            DECred.append(DECg[i])
            V_r.append(V[i])
    
    coords = zip(RAg, DECg)
    coords_blue = zip(RAblue, DECblue)
    coords_red = zip(RAred, DECred)
    
    fits_file = '/home/emilio/MLE/Galaxies/'+gal+'/'+'f.fits'
    hdu = fits.open(fits_file)[0]
    wcs = WCS(hdu.header)
    
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
    
    print 'number of blues:',len(ra_b)
    print 'number of reds:', len(ra_r)
    
    print colored('If the color scale is not ok, change Vmax and Vmin inside the code!', 'yellow')    
    
    Vmax = max(V)-300
    Vmin = min(V)+300
    
#    fig, (ax, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True, subplot_kw=dict(projection=wcs))
#    ax.scatter(ra, dec, c=V, vmin=Vmin, vmax=Vmax)
#    ax.imshow(hdu.data,cmap='gray_r', label='All GCs')
#    ax.set_title('All GCs')    
#    ax2.scatter(ra_b, dec_b, c=V_b, vmin=Vmin, vmax=Vmax)    
#    ax2.imshow(hdu.data,cmap='gray_r', label='Blue GCs')
#    ax2.set_title('Blue GCs')    
#    ax3.scatter(ra_r, dec_r, c=V_r, vmin=Vmin, vmax=Vmax)
#    ax3.imshow(hdu.data,cmap='gray_r', label='Red GCs')
#    ax3.set_title('Red GCs')
#    
#    ax.set_ylabel('RA')
#    ax2.set_xlabel('DEC')
    
    fig = plt.subplot(1,3,1, projection=wcs)
    plt.scatter(ra, dec, c=V, vmin=Vmin, vmax=Vmax)
    plt.imshow(hdu.data,cmap='gray_r', label='All GCs')
    plt.title('All GCs')
    
    fig2= plt.subplot(1,3,2, projection=wcs)
    plt.scatter(ra_b, dec_b, c=V_b, vmin=Vmin, vmax=Vmax)
    cb = plt.colorbar(fraction=0.05, orientation='horizontal')
    cb.set_label('Velocity $(km/s)$')
    plt.imshow(hdu.data,cmap='gray_r', label='Blue GCs')
    plt.title('Blue GCs')
    
    fig3 = plt.subplot(1,3,3, projection=wcs)
    plt.scatter(ra_r, dec_r, c=V_r, vmin=Vmin, vmax=Vmax)
    plt.imshow(hdu.data,cmap='gray_r', label='Red GCs')
    plt.title('Red GCs')
    
    fig.set_ylabel('DEC')
    fig2.set_xlabel('RA')
    
    plt.subplots_adjust(wspace=.0010)
    plt.show()
    
op = raw_input('plot Voronoi diagrams colored by velocity? y/n: ')

if op == 'y':
    plotVor()
    print colored('Plot saved in /plots!', 'yellow')
