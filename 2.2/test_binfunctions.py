import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from bin_functions import bin_areafix, density, binning, plotgcbins, plotgcbins_N, color_dist
import os

# Adjust file names to test the code

print("Density Plot Plotter by Emilio Zanatta\n")
print("Make sure all inputs are in the folder ~/MLE/Galaxies/XXXX, where XXXX is the Galaxie's number.\n")
print("Check also if the catalogs are in the right format:\n")
print("# RA DEC V redder bluer\n")
print("RA and DEC must be in decimal degrees.\n")

op = raw_input('Type GC or PNE: ')
gal = raw_input('Type Galaxy number:')
nbin = raw_input('Enter number of bins: ')
nbin = int(nbin)
cho = raw_input('Enter \'c\' for circular bins or \'e\' for elliptical bins:  ')
bintype = raw_input('Enter A for fixed binsize or N for fixed number of objects per bin: ')

galinput = 'N'+gal+'input.dat'
galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput

if op == 'GC':
    galcat = 'N'+gal+'GC.dat'    
else:
    galcat = 'N'+gal+'PNE.dat'

galcat = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcat

with open(galinput, 'r') as f:
    inp = [x.split(' ')[0] for x in open(galinput).readlines()]
    
c_sep = float(inp[4])

if c_sep == -1000:
    op3 = op
    op = 'UNIMODAL'

RAgal = inp[0]
DECgal = inp[1]
ari = float(inp[2])*u.deg
pa = float(inp[3])*u.deg

galcenter = SkyCoord(RAgal, DECgal)


if op == 'GC':

    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    RA = RA*u.deg
    DEC = DEC*u.deg
    
    V = np.genfromtxt(galcat, usecols=(3,), missing_values=9999., filling_values=0.0)          # velocities
    
    col1 = np.genfromtxt(galcat, usecols=(4,), missing_values=9999., filling_values=0.0)
    col2 = np.genfromtxt(galcat, usecols=(5,), missing_values=9999., filling_values=0.0)
    colour = col1 - col2


    ####### color separating the GCs #######
    n = len(RA)
    
    RAred = []
    DECred = []
    Vred = []
    RAblue = []
    DECblue = []
    Vblue = []
        
                
    
    for i in range(0, n):
        if colour[i] > c_sep:
            RAred.append(RA[i])
            DECred.append(DEC[i])
            Vred.append(V[i])
        else:
            RAblue.append(RA[i])
            DECblue.append(DEC[i])
            Vblue.append(V[i])

    #binning and density        
                    
    
    
    if cho=='c':
        r, rgal, V = bin_areafix(RA, DEC, V, galcenter, pa, ari, nbin, bintype)
        rred, rgal_red, Vred = bin_areafix(RAred, DECred, Vred, galcenter, pa, ari, nbin, bintype)
        rblue, rgal_blue, Vblue = bin_areafix(RAblue, DECblue, Vblue, galcenter, pa, ari, nbin, bintype)
    else:
        r, rgal, V = binning(RA, DEC, V, galcenter, pa, ari, nbin, bintype)
        rred, rgal_red, Vred = binning(RAred, DECred, Vred, galcenter, pa, ari, nbin,bintype)
        rblue, rgal_blue, Vblue = binning(RAblue, DECblue, Vblue, galcenter, pa, ari, nbin,bintype)
        
           
    dens, denslog, poi_err, area, NGC, median, binsize, erro  = density(nbin, r, rgal)
    densred, denslogred, poi_errred, areared, NGCreds, medianreds, binsizereds, erro2 = density(nbin, rred, rgal_red)
    densblue, denslogblue, poi_errblue, areablue, NGCblues, medianblues, binsizeblues, erro3 = density(nbin, rblue, rgal_blue)        
    
    ## outputs for model-elli...py
    
    with open('test_GCdensity-circularbins.dat', 'w') as GC:        #write out the densities
        print >>GC, '# nbin   dens denslog err area NGC median rgal binsize'
        for h in range(0, nbin):
            print >>GC, h, dens[h], denslog[h], poi_err[h], area[h], NGC[h], median[h], rgal[h+1], binsize[h], erro[h]
    with open('test_GCdensity-REDS.dat', 'w') as GC:        #write out the densities
        print >>GC, '# nbin   dens denslog err area NGC median rgal binsize'
        for h in range(0, nbin):
            print >>GC, h, densred[h], denslogred[h], poi_errred[h], areared[h], NGCreds[h], medianreds[h], rgal_red[h+1], binsizereds[h]
    with open('test_GCdensity-BLUES.dat', 'w') as GC:        #write out the densities
        print >>GC, '# nbin   dens denslog err area NGC median rgal binsize' 
        for h in range(0, nbin):
            print >>GC, h, densblue[h], denslogblue[h], poi_errblue[h], areablue[h], NGCblues[h], medianblues[h], rgal_blue[h+1], binsizeblues[h]
    
    ##############################  plots #########################################
    

#    plothist = raw_input('Plot color histogram to check color separation? y/n ')
#    if plothist=='y':    
#        #Colour distribution
#        colour = np.ma.masked_equal(colour, 0.0)
#        hist = plt.subplot(111)
#        plt.hist(colour, histtype='stepfilled', bins=20, color='gray')
#        plt.plot((c_sep, c_sep),(0.0, len(colour)/2), 'k--')
#        #hist.set_xlim([0.1, 1.2])
#        #hist.set_ylim([0.0, 25.0])
#        plt.xlabel('color')
#        plt.ylabel('Number of GC')
#        plt.show()


                  
else:
    #binning and density 

    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    RA = RA*u.deg
    DEC = DEC*u.deg
        
    V = np.loadtxt(galcat, usecols=(3,))
                    
    
    if cho=='A':
        r, rgal, V = bin_areafix(RA, DEC, V, galcenter, pa, ari, nbin, bintype)
    else:
        r, rgal, V = binning(RA, DEC, V, galcenter, pa, ari, nbin, bintype)
    
    
    dens, denslog, poi_err, area, NGC, median, binsize, erro  = density(nbin, r, rgal)
    
    ## outputs for model-elli...py
    
    if (c_sep == -1000 and op3 == 'GC'):
       with open('test_GCdensity-circularbins.dat', 'w') as GC:      #write out the densities
        print >>GC, '# nbin   dens denslog err area NGC median rgal binsize'
        for h in range(0, nbin):
            print >>GC, h, dens[h], denslog[h], poi_err[h], area[h], NGC[h], median[h], rgal[h+1], binsize[h], erro[h]            
    else:
        with open('PNEdensity-circularbins.dat', 'w') as PNE:        #write out the densities
            print >>PNE, '# nbin   dens denslog err area NGC median rgal binsize'
            for h in range(0, nbin):
                print >>PNE, h, dens[h], denslog[h], poi_err[h], area[h], NGC[h], median[h], rgal[h+1], binsize[h], erro[h] 
                
            
with open('current_galaxy.dat', 'w+') as current:     
    print >>current, gal       
########################### PNE version ########################################################                  
                                          
                  
#follow-up to what really matters:       
                  
k = raw_input('plot bins? (only for GC) y/n: ')
if k == 'y':
    if cho=='c':
        plotgcbins()
    if cho=='e':
        plotgcbins_N(ari.value)  

g = raw_input('plot colour vs radius? y/n: ')

if g=='y':
   color_dist(colour)                   
                  
s = raw_input('plot densities? y/n: ')
if s == 'y':
    os.system('python model-ellipse-PNes-plots.py')

