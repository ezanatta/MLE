import numpy as np
import pylab
from astropy import units as u
from astropy.coordinates import SkyCoord
import uncertainties.unumpy as unp

ari = 73*u.deg
ari_rad = ari.to(u.rad)
ari_rad = ari_rad.value
gal_srs_ai = 0.29         #
gal_exp1_ai = 0.66        #  galaxy axis ratio estimated by Galfit 
pa = -86.25*u.deg         #  position angle (from GATOR) 
pa_rad = pa.to(u.rad)
pa_rad = pa_rad.value
galcenter = SkyCoord('09h11m37.5s', '+60d02m14s')
Vsys = 1353 #km/s


with open('2768cat.dat', 'r') as f:
    RA = [line.split()[0] for line in f]
with open('2768cat.dat', 'r') as f:    
    DEC = [line.split()[1] for line in f]

V = np.loadtxt('2768raw_catalog.dat', usecols=(8,))          # velocities
#dV = np.loadtxt('2768raw_catalog.dat', usecols=(4,))        # velocity errors


n = len(RA)         #number of PNe
print 'number of PNE:', n
pne = np.linspace(0,0,n)
ra = np.linspace(0,0,n)
dec = np.linspace(0,0,n)

for i in range(0,n):
    aux = SkyCoord(RA[i], DEC[i])             #reading RA and DEC from PNe
    ra[i] = aux.ra.arcsec
    dec[i] = aux.dec.arcsec
    
        
y_rad = galcenter.dec.radian   # galaxy center DEC coordinate in radians
        
xmc = (ra - galcenter.ra.arcsec)
xm = xmc*(np.cos(y_rad))                       # setting the x coordinates of PNe from its RA value relative to galactic center      
ym = (dec - galcenter.dec.arcsec)             # setting the y coordinates of the PNe from its DEC value relative to galactic center
                     
             
xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # rotating coordinates 
ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    

cos_i = np.cos(ari_rad)               # cosine of inclination angle, used to shrink coordinates
sen_i = np.sin(ari_rad)

xsi = xs*cos_i
ysi = ys/cos_i

r=np.sqrt(((xs**2)*cos_i)+((ys**2)/cos_i))   #distributing the PNe along the radius
cos_phi = xsi/r                              #azimuthal angle cosine & sine
sen_phi = ysi/r  

rvel = r

##############  binning section begins   ######################################################################

t = 0
n = r.size   #lenght of r --- number of objects to iterate on the following loops
for i in range(0, n):
    for j in range(i+1,n):
        if r[i] > r[j]:
               t = r[i]
               r[i] = r[j]
               r[j] = t
                             
nbin = 4                     #number of bins
rgal = np.linspace(0,0,nbin+1)       #rgal contains the limits of each bin

for h in range(1, nbin):
    rgal[h] = r[int(n*h/nbin)]     #binning
rgal[nbin] = r[n-1]
rgal[0] = r[0]


NPNe = np.zeros(nbin)    #number of PNe in each bin
rmed = [[] for i in range(0,nbin)]
Vbin = [[] for i in range(0,nbin)]
for h in range(0, nbin):                              
     for i in range(0, n):
         if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
             NPNe[h] = NPNe[h] + 1
             rmed[h].append(r[i])
             
#calculation of mean streaming motion (V_r):

V_r = np.linspace(0,0,n)
for i in range(0,n):
    V_r[i] = (V[i]-Vsys)/(np.sin(ari)*cos_phi[i])
    print cos_phi[i], np.sin(ari)
             
for h in range(0,nbin):
    for i in range(0,n):
        if (rgal[h+1] >= rvel[i] and rvel[i] > rgal[h]): 
            Vbin[h].append(V_r[i])
            
binsize = np.linspace(0,0,nbin)    
for h in range(0,nbin):
    binsize[h] = rgal[h+1]-rgal[h]  
binsize = binsize/2
            
Vavr_bin = np.linspace(0,0,nbin)            
for h in range(0,nbin):
    Vsum_bin = sum(Vbin[h])
    Vavr_bin[h] = Vsum_bin/NPNe[h]
        
#median value of the bins
radius = [[] for i in range(0,nbin)]
for i in range(0,nbin):
    radius[i] = np.linspace(rgal[i], rgal[i+1], 1000)   
median = np.linspace(0,0,nbin)
medianPNE = np.linspace(0,0,nbin)
for i in range(0,nbin):
    median[i] = np.median(radius[i])
    medianPNE[i] = np.median(rmed[i])
    
with open('velocity_plots_results.dat', 'w') as f:
    for h in range(0,nbin):
        print >>f, h, NPNe[h], Vavr_bin[h]
    
pylab.errorbar(medianPNE, Vavr_bin, xerr=binsize, marker = 'o', color = 'red', linestyle='none')
pylab.show()
