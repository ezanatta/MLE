import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import binning
import os.path

class pf(float):
    def __repr__(self):
        return "%0.2f" % self
        
font = {'family': ['Computer Modern'],
        'color':  'black',
        'weight': 'normal',
        'size': 20,
        }        

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
    
def read_catalog_stars(starcat):
    with open(starcat, 'r') as f:
        RAs = [line.split()[0] for line in f]               
    with open(starcat, 'r') as f:    
        DECs = [line.split()[1] for line in f]
    
    for i in range(0, len(RAs)):
        c = SkyCoord(RAs[i], DECs[i])
        RAs[i], DECs[i] = c.ra.degree, c.dec.degree
 
    RAs = np.asarray(RAs)
    DECs = np.asarray(DECs)

    Vs = np.loadtxt(starcat, usecols=(2,))
    
    return RAs, DECs, Vs

    
def pne_radius(RA, DEC, gc, pa, i):

     pa = pa
     pa_rad = pa.to(u.rad)
     pa_rad = pa_rad.value  
     
     ra = np.linspace(0,0,len(RA))
     dec = np.linspace(0,0,len(RA))


     cos_i = np.cos(i.to(u.rad))
     yrad = np.cos(gc.dec.radian)

    # print 'cos(DECg) = ', yrad
     
     for i in range(0,len(RA)):
         aux = SkyCoord(RA[i]*u.deg, DEC[i]*u.deg)             #reading RA and DEC from GC
         ra[i] = aux.ra.arcsec
         dec[i] = aux.dec.arcsec  
         
     xm = (ra - gc.ra.arcsec)
     ym = (dec - gc.dec.arcsec)
     xm = xm*yrad
     
     xsi = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)
     ysi = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)
     
     xsi = xsi*np.sqrt(cos_i)
     ysi = ysi/np.sqrt(cos_i)
     
     return -xsi, ysi
     
def atokpc(R, d):
    from astropy.cosmology import WMAP9 as cosmo
    if gal == '2768':
        z = 0.004513
    elif gal == '7457':
        z = 0.002815
    elif gal == '3115':
        z = 0.002212
    
    Dkpc = cosmo.kpc_proper_per_arcmin(z).value*(R/60)
    
    return Dkpc    
     

gal = raw_input('Type Galaxy number:')
#op = raw_input('Do we have GC? y/n: ')
#op2 = raw_input('Do we have PNe? y/n: ')
op='y'
op2='y'

galinput = 'N'+gal+'input.dat'
galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput

if op == 'y' and op2 == 'y':
    galcatgc = 'N'+gal+'GC.dat'    
    galcatpne = 'N'+gal+'PNE.dat'
    #galcatstars = 'N'+gal+'STARSProctor08.dat'
    
    galcatgc = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatgc
    galcatpne = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatpne
    #galcatstars = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatstars
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]
        
    c_sep = float(inp[4])
    
    RAgal = inp[0]
    DECgal = inp[1]
    i = float(inp[2])*u.deg
    pa = float(inp[3])*u.deg
    Vsys = float(inp[10])
    d = float(inp[11])*u.mpc
    
    galcenter = SkyCoord(RAgal, DECgal)
    
    RA, DEC, V, col = read_catalog_gc(galcatgc)
    RAp, DECp, Vp = read_catalog_pne(galcatpne) 
    
    Vp = Vp-Vsys
    
    nbin = 6
    nbinp = 8
    
    rp, rbinp, Vordp = binning(RAp, DECp, Vp, galcenter, pa, i, nbinp, 'N')
    r, rbin, Vord = binning(RA, DEC, V, galcenter, pa, i, nbin, 'N')
    
    
    RA, DEC, V, col = read_catalog_gc(galcatgc)
    RAp, DECp, Vp = read_catalog_pne(galcatpne) 
    Vp = Vp-Vsys
    
    op = raw_input('which backup? ')
    
    if os.path.exists('../2.5/N'+gal+'/backup'+op+'/prob_out.dat'):
        
        RA = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(0,))
        DEC = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(1,))
        fb = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(5,))
        vfb = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(4,))
        #xs = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(2,))
        #ys = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(3,))
        colf = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(7,))
        like = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/prob_out.dat', usecols=(8,))
        
    elif os.path.exists('../2.5/N'+gal+'/backup'+op+'/unliky.dat'):  
        
        RA = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(0,))
        DEC = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(1,))
        fb = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(6,))
        vfb = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(5,))
        #xs = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(2,))
        #ys = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(3,))
        colf = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(9,))
        like = np.loadtxt('../2.5/N'+gal+'/backup'+op+'/unliky.dat', usecols=(13,))

    opp = 0.0*u.deg
    while opp!='-1':    
    
        px, py = pne_radius(RAp.value, DECp.value, galcenter, pa, i)
        gx, gy = pne_radius(RA, DEC, galcenter, pa, i)            
        
        #plt.scatter(gx, gy, c=V, cmap='rainbow', vmin=np.min(V), vmax=np.max(V))
        #plt.show()
        
        print 'current PA: ', pa
        #opp = raw_input('try another PA? Enter value or type ''-1'' to continue: ')
        opp = '-1'
        pa = np.float(opp)*u.deg        
#    ys = ys/np.sqrt(np.cos(i.to(u.rad)))

#    gx = map(pf, gx)
#    gy = map(pf, gy)
#
    gx = atokpc(gx, d)
    gy = atokpc(gy, d)
    px = atokpc(px, d)
    py = atokpc(py, d)
    
    xs = gx
    ys = gy

    vfb = vfb-Vsys
    
    xrej = []
    yrej = []
    colrej = []
    vrej = []
        
    for i in range(0, len(like)):
        if like[i] == 5.00:
            xrej.append(xs[i])    
            yrej.append(ys[i])    
            colrej.append(colf[i])
            vrej.append(vfb[i])     

    xrej = np.asarray(xrej)
    yrej = np.asarray(yrej)
    colrej = np.asarray(colrej)
    vrej = np.asarray(vrej)                
    

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
            xsred.append(ys[i])
        else:
            rblue.append(np.sqrt(xs[i]**2+ys[i]**2))
            Vblue.append(vfb[i])  
            fblue.append(fb[i])
            xsblue.append(ys[i])
            
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

    print(len(rdiskblue)+len(rdiskred))
    print(len(rsphblue)+len(rsphred))
    
    if gal == '2768':
        xmaxi = 60
        colmax = 1.2
        colmin = 0.2
        col_label = '(Rc-z)'
    elif gal == '3115':
        xmaxi = 45
        col_label='(g-i)'
        colmin=0.6
        colmax=1.4
    elif gal == '7457':
        xmaxi = 17

    if c_sep != -1000:
        #plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
#        plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
#        plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
        plt.plot(np.sqrt(px**2+py**2), Vp, marker='^', markersize=6, markerfacecolor='lightgreen', markeredgecolor='lightgreen', linestyle='None', label= 'PNe')
        plt.plot(rdiskred, vdiskred, markersize=8,marker='D', markerfacecolor='None', markeredgecolor='red', linestyle='None')
        plt.plot(rdiskblue, vdiskblue, markersize=8, marker='D', markerfacecolor='None', markeredgecolor='blue', linestyle='None')
        plt.plot(rsphred, vsphred, markersize=8, marker='o', color='red', linestyle='None')
        plt.plot(rsphblue, vsphblue, markersize=8, marker='o', color='blue', linestyle='None') 
        plt.plot(np.sqrt(xrej**2+yrej**2), vrej, markersize=10, color='black', marker='h',linestyle='None', label='Rejected GC')
        plt.xlim(0, xmaxi)
        plt.ylabel('$\Delta$V $(km/s)$', fontdict=font)
        plt.xlabel('R $(kpc)$', fontdict=font)
        plt.plot((0, xmax),(0, 0), ':', label='Galaxy Systemic Velocity', color='black')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.savefig(gal+'phasespace.png', dpi=300)
        plt.show()  
        
        #    plt.scatter(xs, vfb, c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0)
        plt.plot(py, Vp, marker='o', markersize=6, markerfacecolor='lightgreen', markeredgecolor='lightgreen', linestyle='None', label= 'PNe')
        #plt.plot(sx, Vs, marker='o', markersize=4, markerfacecolor='lightgreen', markeredgecolor='lightgreen', linestyle='None')
        plt.plot(xdiskred, vdiskred, markersize=8, marker='D', markerfacecolor='None', markeredgecolor='red', linestyle='None')
        plt.plot(xdiskblue, vdiskblue, markersize=8, marker='D', markerfacecolor='None', markeredgecolor='blue', linestyle='None')
        plt.plot(xsphred, vsphred, markersize=8, marker='o', color='red', linestyle='None')
        plt.plot(xsphblue, vsphblue, markersize=8, marker='o', color='blue', linestyle='None')      
        plt.plot((-700, 700),(0, 0), ':', label='Galaxy Systemic Velocity', color='black')
        plt.plot((0, 0),(-700, 700), ':', color='black')
        plt.xlim(-xmaxi, xmaxi)
        plt.ylim(-500, 500)
        plt.plot(yrej, vrej, color='black', markersize=10,marker='h',linestyle='None', label='Rejected GC') 
        plt.ylabel('$\Delta$V $(km/s)$', fontdict=font)
        plt.xlabel('$\Delta$x $(kpc)$', fontdict=font)
    #    cb = plt.colorbar(fraction=0.05)
    #    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.savefig(gal+'crossspace.png', dpi=300)
        plt.show()
        
        fig3 = plt.subplot(1, 1, 1)
        plt.scatter(colf, np.sqrt(xs**2+ys**2), c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0, s=50)
        for i in range(0, len(colrej)):
            fig3.plot(colrej[i], np.sqrt(xrej[i]**2+yrej[i]**2), color='red', marker='x', markersize=10,linestyle='None',mew=3, label='Rejected GC')    
        fig3.set_xlim(colmin, colmax)   
        fig3.set_ylim(0, xmaxi)    
        fig3.set_xlabel(col_label+'(mag)', fontsize=20)
        fig3.set_ylabel('R$(kpc)$', fontsize=20)
        cb = plt.colorbar(fraction=0.05)
        cb.set_label('$L_{Sph}(v_{i}, f_{i})$', fontsize=20)
        plt.savefig(gal+'rcolf.png', dpi=300)
        plt.show()
    
    else:
        #plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
        plt.plot(np.sqrt(px**2+py**2), Vp, marker='^', markersize=6, markerfacecolor='lightgreen', markeredgecolor='lightgreen', linestyle='None', label= 'PNe')
        plt.plot(rdiskred, vdiskred, marker='o', markersize=8,markerfacecolor='None', markeredgecolor='purple', linestyle='None')
        plt.plot(rsphred, vsphred, marker='o', markersize=8,color='purple', linestyle='None')
        #plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
        #plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
    #    plt.plot(rp, Vordp, marker='o', markerfacecolor='None', markeredgecolor='green', linestyle='None', label= 'PNe')
        plt.ylabel('$\Delta$V $(km/s)$', fontdict=font)
        plt.xlabel('R $(kpc)$', fontdict=font)
        plt.xlim(0, xmaxi)
        plt.plot(np.sqrt(xrej**2+yrej**2), vrej, markersize=10, color='black', marker='h',linestyle='None', label='Rejected GC')
        plt.plot((0, xmax),(0, 0), '--', label='Galaxy Systemic Velocity', color='black')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.savefig(gal+'phasespace.png', dpi=300)
        plt.show()  
          
    
    #    plt.scatter(xs, vfb, c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0)
        plt.plot(py, Vp, marker='^', markersize=6, markerfacecolor='lightgreen', markeredgecolor='lightgreen', linestyle='None', label= 'PNe')
        plt.plot(np.asarray(xdiskred), vdiskred, markersize=8,marker='o', markerfacecolor='None', markeredgecolor='purple', linestyle='None')
#        plt.plot(xdiskblue, vdiskblue, marker='o', markerfacecolor='None', markeredgecolor='blue', linestyle='None')
        plt.plot(np.asarray(xsphred), vsphred, markersize=8,marker='o', color='purple', linestyle='None')
#        plt.plot(xsphblue, vsphblue, marker='o', color='blue', linestyle='None')      
        plt.plot((-700, 700),(0, 0), ':', label='Galaxy Systemic Velocity', color='black')
        plt.plot((0, 0),(-700, 700), ':', color='black')
        plt.xlim(-8, 8)
        plt.ylim(-150, 200)
        plt.plot(yrej, vrej, color='black', markersize=10,marker='x',linestyle='None', label='Rejected GC') 
        plt.ylabel('$\Delta$V $(km/s)$', fontdict=font)
        plt.xlabel('$\Delta$x $(kpc)$', fontdict=font)
    #    cb = plt.colorbar(fraction=0.05)
    #    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
        #plt.legend(loc='upper right', prop={'size':10})
        plt.savefig(gal+'crossspace.png', dpi=300)
        plt.show()
        
        fig3 = plt.subplot(1, 1, 1)
        plt.scatter(colf, np.sqrt(xs**2+ys**2), c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0, s=50)
        for i in range(0, len(colrej)):
            fig3.plot(colrej[i], np.sqrt(xrej[i]**2+yrej[i]**2),color='red', marker='x', markersize=10,linestyle='None',mew=4, label='Rejected GC')    
        fig3.set_xlabel('($R_{c}$-z)(mag)')
        fig3.set_ylabel('R$(kpc)$')
        cb = plt.colorbar(fraction=0.05)
        cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
        plt.savefig(gal+'rcolf.png', dpi=300)
        plt.show()
        