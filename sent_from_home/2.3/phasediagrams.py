import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import binning

def read_catalog_gc(galcat):
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    RA = RA*u.deg
    DEC = DEC*u.deg    
    
    V = np.loadtxt(galcat, usecols=(3,))        
    
    col1 = np.loadtxt(galcat, usecols=(4,))
    col2 = np.loadtxt(galcat, usecols=(5,))
    colour = col1 - col2
    
    return RA, DEC, V, colour  
    
def read_catalog_pne(pnecat):
    RAp = np.loadtxt(pnecat, usecols=(1,))
    DECp = np.loadtxt(pnecat, usecols=(2,))
    RAp = RAp*u.deg
    DECp = DECp*u.deg
    
    Vp = np.loadtxt(pnecat, usecols=(3,))
    
    return RAp, DECp, Vp

gal = raw_input('Type Galaxy number:')
op = raw_input('Do we have GC? y/n: ')
op2 = raw_input('Do we have PNe? y/n: ')

galinput = 'N'+gal+'input.dat'
galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput

if op == 'y' and op2 == 'y':
    galcatgc = 'N'+gal+'GC.dat'    
    galcatpne = 'N'+gal+'PNE.dat'
    
    galcatgc = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatgc
    galcatpne = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatpne
    
    fb = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(6,))
    vfb = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(5,))
    xs = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(2,))
    ys = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(3,))
    colf = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(9,))
    like = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/unliky.dat', usecols=(10,))
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]
        
    c_sep = float(inp[4])
    
    RAgal = inp[0]
    DECgal = inp[1]
    i = float(inp[2])*u.deg
    pa = float(inp[3])*u.deg
    Vsys = float(inp[10])
    
    galcenter = SkyCoord(RAgal, DECgal)
    
    RA, DEC, V, col = read_catalog_gc(galcatgc)
    RAp, DECp, Vp = read_catalog_pne(galcatpne) 
    
    vfb = vfb-Vsys
    
    nbin = 6
    nbinp = 8
    
    rp, rbinp, Vordp = binning(RAp, DECp, Vp, galcenter, pa, i, nbinp, 'N')
    r, rbin, Vord = binning(RA, DEC, V, galcenter, pa, i, nbin, 'N')
    
    rred = []
    rblue = []
    Vred = []
    Vblue = []    
    lsph = []
    fred = []
    fblue = []
    xsred = []
    xsblue = []
    
    for i in range(0, len(r)):
        if colf[i] >= c_sep:
            rred.append(np.sqrt(xs[i]**2+ys[i]**2))
            Vred.append(vfb[i])
            fred.append(fb[i])
            xsred.append(xs[i])
        else:
            rblue.append(np.sqrt(xs[i]**2+ys[i]**2))
            Vblue.append(vfb[i])  
            fblue.append(fb[i])
            xsblue.append(xs[i])
            
    #selecting objects within the colors in disk and spheroid subpop.
    rsphred = []
    vsphred = []
    rsphblue = []
    vsphblue = []
    rdiskred = []
    rdiskblue = []
    vdiskblue = []
    vdiskred = []    

            
    xsphred = []
    xsphblue = []
    xdiskred = []
    xdiskblue = []           
            
            
    for i in range(0, len(rred)):
       if fred[i] > 0.5:
           rsphred.append(rred[i])
           vsphred.append(Vred[i])
           xsphred.append(xsred[i])
       else:
           rdiskred.append(rred[i])
           vdiskred.append(Vred[i])
           xdiskred.append(xsred[i])
    for i in range(0, len(rblue)):
       if fblue[i] > 0.5:
           rsphblue.append(rblue[i])
           vsphblue.append(Vblue[i])
           xsphblue.append(xsblue[i])
       else:
           rdiskblue.append(rblue[i])
           vdiskblue.append(Vblue[i])           
           xdiskblue.append(xsblue[i])
            
    xmax = max(r)+200
    
    if c_sep != -1000:
        #plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
#        plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
#        plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
        plt.plot(rdiskred, vdiskred, marker='D', markerfacecolor='None', markeredgecolor='red', linestyle='None')
        plt.plot(rdiskblue, vdiskblue, marker='D', markerfacecolor='None', markeredgecolor='blue', linestyle='None')
        plt.plot(rsphred, vsphred, marker='o', color='red', linestyle='None')
        plt.plot(rsphblue, vsphblue, marker='o', color='blue', linestyle='None')        
        
    #    plt.plot(rp, Vordp, marker='o', markerfacecolor='None', markeredgecolor='green', linestyle='None', label= 'PNe')
        plt.ylabel('V $(km/s)$')
        plt.xlabel('R $(arcsec)$')
        plt.plot((0, xmax),(0, 0), '--', label='Galaxy Systemic Velocity', color='black')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.show()   
    else:
        #plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
        plt.plot(rdiskred, vdiskred, marker='D', markerfacecolor='None', markeredgecolor='purple', linestyle='None')
        plt.plot(rsphred, vsphred, marker='o', color='purple', linestyle='None')
        #plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
        #plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
    #    plt.plot(rp, Vordp, marker='o', markerfacecolor='None', markeredgecolor='green', linestyle='None', label= 'PNe')
        plt.ylabel('V $(km/s)$')
        plt.xlabel('R $(arcsec)$')
        plt.plot((0, xmax),(0, 0), '--', label='Galaxy Systemic Velocity', color='black')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.show()  
    
    plt.scatter(np.sqrt(xs**2+ys**2), vfb, c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0)
    #plt.plot(rp, Vordp, marker='o', markerfacecolor='None', markeredgecolor='green', linestyle='None', label= 'PNe')
    plt.ylabel('$\Delta$V $(km/s)$')
    plt.xlabel('R $(arcsec)$')
    plt.plot((0, xmax),(0, 0), '--', label='Galaxy Systemic Velocity')
    cb = plt.colorbar(fraction=0.05)
    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
    #plt.legend(loc='upper right', prop={'size':10})
    plt.show()
    
    xrej = []
    yrej = []
    colrej = []
    
    for i in range(0, len(like)):
        if np.abs(like[i]) > 10.00:
            xrej.append(xs[i])      
            yrej.append(ys[i])    
            colrej.append(vfb[i])         
    
#    plt.scatter(xs, vfb, c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0)
    
    plt.plot(xdiskred, vdiskred, marker='D', markerfacecolor='None', markeredgecolor='red', linestyle='None')
    plt.plot(xdiskblue, vdiskblue, marker='D', markerfacecolor='None', markeredgecolor='blue', linestyle='None')
    plt.plot(xsphred, vsphred, marker='o', color='red', linestyle='None')
    plt.plot(xsphblue, vsphblue, marker='o', color='blue', linestyle='None')      
    plt.plot((-400, 400),(0, 0), '--', label='Galaxy Systemic Velocity', color='blue')
    plt.plot((0, 0),(min(vfb), max(vfb)), '--', color='cyan')
    #plt.plot(xrej, colrej, 'ro',linestyle='None', label='Rejected GC') 
    plt.ylabel('$\Delta$V $(km/s)$')
    plt.xlabel('$\Delta$x $(arcsec)$')
#    cb = plt.colorbar(fraction=0.05)
#    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
    #plt.legend(loc='upper right', prop={'size':10})
    plt.show()

if op == 'y' and op2 != 'y':
    galcatgc = 'N'+gal+'GC.dat'    
    
    galcatgc = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatgc
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]
        
    c_sep = float(inp[4])
    
    RAgal = inp[0]
    DECgal = inp[1]
    i = float(inp[2])*u.deg
    pa = float(inp[3])*u.deg
    Vsys = float(inp[10])
    
    galcenter = SkyCoord(RAgal, DECgal)
    
    RA, DEC, V, col = read_catalog_gc(galcatgc)
    
    nbin = 6

    r, rbin, Vord = binning(RA, DEC, V, galcenter, pa, i, nbin, 'N')
    
    xmax = max(r)+200
    
    plt.plot(r, Vord, marker='o', markeredgecolor='purple', markerfacecolor='None', linestyle='None', label='GC')
    plt.ylabel('V $(km/s)$')
    plt.xlabel('R $(arcsec)$')
    plt.plot((0, xmax),(Vsys, Vsys), '--', label='Galaxy Systemic Velocity')
    plt.legend(loc='upper right', prop={'size':10})
    plt.show()

if op != 'y' and op2 == 'y':   
    galcatpne = 'N'+gal+'PNE.dat'
    
    galcatpne = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatpne
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]
        
    c_sep = float(inp[4])
    
    RAgal = inp[0]
    DECgal = inp[1]
    i = float(inp[2])*u.deg
    pa = float(inp[3])*u.deg
    Vsys = float(inp[10])
    
    galcenter = SkyCoord(RAgal, DECgal)
    
    RAp, DECp, Vp = read_catalog_pne(galcatpne) 
    
    nbinp = 8
    
    rp, rbinp, Vordp = binning(RAp, DECp, Vp, galcenter, pa, i, nbinp, 'N')
    
    xmax = max(rp)+200
    
    plt.plot(rp, Vordp, marker='o', color='lightgreen', linestyle='None', label= 'PNe')
    plt.ylabel('V $(km/s)$')
    plt.xlabel('R $(arcsec)$')
    plt.plot((0, xmax),(Vsys, Vsys), '--', label='Galaxy Systemic Velocity')
    plt.legend(loc='upper right', prop={'size':10})
    plt.show()    



