#import matplotlib
#matplotlib.use('Agg')               #if remotely connect decomment this
import numpy as np
import pylab
from matplotlib.patches import Ellipse
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
                             
nbin = 8                     #number of bins
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
             
#median value of the bins
radius = [[] for i in range(0,nbin)]
for i in range(0,nbin):
    radius[i] = np.linspace(rgal[i], rgal[i+1], 1000)   
median = np.linspace(0,0,nbin)
medianPNE = np.linspace(0,0,nbin)
for i in range(0,nbin):
    median[i] = np.median(radius[i])
    medianPNE[i] = np.median(rmed[i])
    
corr = np.loadtxt('PNeextrapolation.dat', usecols=(1,))   #completness correction
N = np.zeros(nbin)                                      #number of PNe + Completness
for h in range(0,nbin):
    N[h] = NPNe[h] + NPNe[h]*(1-(corr[h]))

area = []
for h in range(0, nbin):
    area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
    
binsize = np.linspace(0,0,nbin)    
for h in range(0,nbin):
    binsize[h] = rgal[h+1]-rgal[h]  
binsize = binsize/2
    
    
dens = np.linspace(0,0,nbin)   #density per bin (N/area)
dens2 = np.linspace(0,0,nbin)
for h in range(0, nbin):
    dens[h] = (N[h]/area[h])
    dens2[h] = NPNe[h]/area[h]
    
#error in density

poi_err = (np.sqrt(N)/area)
denserr = unp.uarray(dens, poi_err)
dens2err = unp.uarray(dens2, poi_err)    

denslog = -2.5*np.log10(dens)
denslog2 = -2.5*np.log10(dens2)

denslogerr =  -2.5*unp.log10(denserr)   #with correction
denslog2err = -2.5*unp.log10(dens2err)  #without correction

poi_err = unp.std_devs(denslogerr)
poi_err2 = unp.std_devs(denslog2err)


    
with open('PNEdensity-circularbins.dat', 'w') as pne:        #write out the densities
    for h in range(0, nbin):
        print >>pne, h, dens2[h], denslog2[h], poi_err2[h], area[h], NPNe[h], median[h], rgal[h+1], binsize[h]

with open('PNElog-circular.dat', 'w') as catalog:    
    for i in range(0,n):
        print >>catalog, xs[i], ys[i]

##############################  plots #########################################
#some test plots to understand better what i'm doing in the end
phi = np.arccos(cos_phi)               # azimuthal angle
 
#color coding the PNe with respect to velocities
pylab.scatter(-xsi, ysi, c=V)
pylab.plot(0,0,'rx', markersize=8, label='galaxy center')
cb = pylab.colorbar()
cb.set_label('Radial Velocity')
pylab.xlabel('x (arcsec)')
pylab.ylabel('y (arcsec)') 
pylab.show()

pylab.plot(-xsi, ysi, 'go')
pylab.xlabel('x (arcsec)')
pylab.ylabel('y (arcmin)')
circles = pylab.subplot(111)
for h in range(0, nbin+1):
    cir = Ellipse((0,0), width=rgal[h], height=rgal[h], fill=False, label='Bins')
    circles.add_artist(cir)                        
pylab.legend(loc='lower left')
pylab.show()     

#velocity distribution
pylab.hist(V, bins=100)
pylab.xlabel('Velocity')
pylab.ylabel('Number of PNe')
pylab.show()


#pylab.errorbar(median, denslog2, yerr=poi_err, xerr=binsize, marker='o', color='blue', linestyle='none', mfc='None' )
pylab.plot(median, denslog, marker='o', color='blue', linestyle='none', mfc='None', label='With Correction')
pylab.plot(median, denslog2, marker='D', color='green', linestyle='none', label='No Correction')
pylab.xlabel('Radius')
pylab.ylabel('PNe Density')
pylab.legend(loc='lower left')
pylab.plt.gca().invert_yaxis()
pylab.show()    
        
pylab.plot(rvel, V, 'ro')
pylab.xlabel('R(arcseconds)')
pylab.ylabel('V (km/s)')
pylab.plt.gca().invert_yaxis()
pylab.show()        
