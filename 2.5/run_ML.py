import sys
import os
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import galaxy_inputs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
from termcolor import colored
import shutil

##############################################################################
# This code is used as way to run the ML.f code in a python sequence. 
# First we need to create the input files for ML.f
# Then we proceed to compile and run ML.f
# TO DO: The outputs are then send to plots_from_like.py to be plotted, if asked.
#
##############################################################################

def read_inputs(galinp):
    inp = [x.split(' ')[0] for x in open(galinp).readlines()]
        
    RAgal = inp[0]
    DECgal = inp[1]
    ari = float(inp[2])
    pa = float(inp[3])
    c_sep = float(inp[4])
    n = float(inp[5])
    re = float(inp[6])
    h = float(inp[7])
    mue = float(inp[8])
    mud = float(inp[9])
    Vsys = float(inp[10])
    d = float(inp[11])    
    
    return  RAgal, DECgal, ari, pa, c_sep, n, re, h, mue, mud, Vsys, d
    
def count_col(cat, c_sep):
    import numpy as np    
    
    col1 = np.loadtxt(cat, usecols=(4,))  
    col2 = np.loadtxt(cat, usecols=(5,))
    col = col1-col2
    
    cont=0
    for item in col:
        if item >= c_sep:
            cont=cont+1
    print 'number of reds:', cont
    print 'number of blues:', len(col)-cont
    print 'number of all: ', len(col)            
    
    
def run_ml():
    os.system('gfortran ML.GC.red.disk.1023.f')
    os.system('./a.out')    
            
    print colored('ML.f finished!', 'yellow')
    
    op2 = raw_input('edit files? y/n ')
    
    if op2 == 'y':
        op3 = raw_input('which color population was run this time? all/red/blue ')
        
        if op3 == 'all':
            os.system('mv likelihood.dat like_.dat')
            os.system('mv bin.dat bin_.dat')
            os.system('mv error.dat error_.dat')
            
        if op3 == 'red':
            os.system('mv likelihood.dat like_red.dat')
            os.system('mv bin.dat bin_red.dat')
            os.system('mv error.dat error_red.dat')
            
        if op3 == 'blue':
            os.system('mv likelihood.dat like_blue.dat')
            os.system('mv bin.dat bin_blue.dat')
            os.system('mv error.dat error_blue.dat')
        
        print colored('finished editing files!', 'red')
        
def run_prob_f(ver):
    nbins = sum(1 for line in open('/home/emilio/MLE/pne_likes/bin_'+gal+'.dat'))
    
    os.system('cp /home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat like_pne.dat')
    os.system('cp /home/emilio/MLE/pne_likes/bin_'+gal+'.dat bin_pne.dat')
    os.system('cp /home/emilio/MLE/Galaxies/'+gal+'/bin_rejpne.dat bin_rejpne.dat')
    
    with open('prob_input.dat', 'w+') as f:
        print >>f, xp, yp, xss, yss, ps, ra, dec, PA, i, V, nbins 
        
    if ver == 'd':
        os.system('gfortran -o prob_a.out prob_bd_3115.f')
        os.system('./prob_a.out')
    if ver =='b':
        os.system('gfortran -o prob_b.out prob_bd_bulge.f')
        os.system('./prob_b.out')
    
    print colored('prob_bulge finished!', 'yellow')
    
def plots_from_prob():
    import matplotlib.pyplot as plt
    import numpy as np    
    
    f = np.loadtxt(f_file, usecols=(0,))
    fb = np.loadtxt('unliky.dat', usecols=(6,))
    xs = np.loadtxt('unliky.dat', usecols=(2,))
    ys = np.loadtxt('unliky.dat', usecols=(3,))
    col = np.loadtxt('unliky.dat', usecols=(9,))
    like = np.loadtxt('unliky.dat', usecols=(13,))
    
    xrej = []
    yrej = []
    colrej = []
    
    print like
    
    for i in range(0, len(like)):
        if np.abs(like[i]) > 0.00:
            xrej.append(xs[i])      
            yrej.append(ys[i])    
            colrej.append(col[i])                    
    
    fig = plt.subplot(1, 1, 1)
    plt.hist(f, bins=20, color='gray')
    fig.set_xlabel('$f_{i}$')
    fig.set_ylabel('N')
    plt.show()
    fig2 = plt.subplot(1, 1, 1)
    plt.hist(fb, bins=20, color='gray')
    fig2.set_xlabel('$L_{Sph}(v_{i}, f_{i})$')
    fig2.set_ylabel('N')
    plt.show()
    fig3 = plt.subplot(1, 1, 1)
    plt.scatter(col, np.sqrt(xs**2+ys**2), c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0, s=50)
    for i in range(0, len(colrej)):
        fig3.plot(colrej[i], np.sqrt(xrej[i]**2+yrej[i]**2), 'ro',linestyle='None', label='Rejected GC')    
    fig3.set_xlabel('(B-V)(mag)')
    fig3.set_ylabel('R$(arcsec)$')
    cb = plt.colorbar(fraction=0.05)
    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
    plt.show()
    
#    plt.hist([f, fb], bins=20, histtype='barstacked', color=['gray', 'white'])
#    plt.xlabel('$f_{i}$')
#    plt.ylabel('N')    
#    plt.show()
   
        
###############'''' main code ''''############################################
    
gal = raw_input('Enter the galaxy number: ')
galinput = galaxy_inputs(gal)
galcat = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'GC.dat' 

with open('current_galaxy.dat', 'w+') as current:     
    print >>current, gal 
    
RA, DEC, i, PA, c, n, Re, Rd, ue, ud, V, d = read_inputs(galinput)
PA = 120

coor = SkyCoord(RA, DEC)

ra = coor.ra.degree
dec = coor.dec.degree

input_name = '/home/emilio/MLE/2.5/input_params.dat'

ps = 1.01  # 2MASS pixel scale 
    
#print colored('Remember to change the center of the fits image inside the code, for every new galaxy!', 'yellow')    

#galcenter = SkyCoord.from_name('NGC'+gal)
#hdulist = fits.open('/home/emilio/MLE/Galaxies/'+gal+'/fits-null.fits')
#hwcs = wcs.WCS(hdulist[0].header)
#center = wcs.utils.skycoord_to_pixel(galcenter, hwcs)
    
xp = 0.0
yp = 0.0
    
xss = xp
yss = yp     #it seems that this variables are actually irrelevant at the moment
    
with open(input_name, 'w+') as f:
    print >>f, xp, yp, xss, yss, ps, ra, dec, PA, i, V, c 
        
gal_cat = '/home/emilio/MLE/Galaxies/'+gal+'/N'+gal+'GC.dat'
f_file = '/home/emilio/MLE/Galaxies/'+gal+'/N'+gal+'_f.GC.dat'
cat_name = '/home/emilio/MLE/2.5/input_catalogue.dat'
f_name = '/home/emilio/MLE/2.5/ffile.dat'
lcut_file = '/home/emilio/MLE/pne_likes/lcut_'+gal+'.dat'
lcut_name = '/home/emilio/MLE/2.5/lcut.dat'

shutil.copyfile(gal_cat, cat_name)
shutil.copyfile(f_file, f_name)
shutil.copyfile(lcut_file, lcut_name)

count_col(gal_cat, c)

op = raw_input('Configuration files ready! proceed to ML.f? y/n ')

if op == 'y':
    run_ml()
print ('\n')
print('\n')
op5 = raw_input('run ML.f again for another population? y/n ')
if op5 == 'y':
    run_ml()
    op6 = raw_input('run ML.f again for yet another population? y/n ')
    if op6 == 'y':
        run_ml()   
        print colored('Nothing else to be done with ML.f\n', 'yellow')
        
op4 = raw_input('Proceed to plots_from_like? y/n ')

if op4 == 'y':
    os.system('python plots_from_runML.py')
else:
    print colored('bye!', 'white')
    
op5 = raw_input('Proceed to prob_bd_bulge.f? y/n ')

if op5 == 'y':
    ver = raw_input('Enter d for disk+bulge or b for only bulge: ')
    run_prob_f(ver)
    plots_from_prob()
