#import matplotlib
#matplotlib.use('Agg')               #if remotely connect decomment this
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
pa = -86.25*u.deg         #position angle (from GATOR) 
pa_rad = pa.to(u.rad)
pa_rad = pa_rad.value
galcenter = SkyCoord('09h11m37.5s', '+60d02m14s')


#with open('2768GC-complete_clean_dec.dat', 'r') as f:
    #RA = [line.split()[0] for line in f]
#with open('2768GC-complete_clean_dec.dat', 'r') as f:    
    #DEC = [line.split()[1] for line in f]

RA = np.loadtxt('2768GC-complete_clean_dec.dat', usecols=(1,))
DEC = np.loadtxt('2768GC-complete_clean_dec.dat', usecols=(2,))

RA = RA*u.deg
DEC = DEC*u.deg

V = np.loadtxt('2768GCcatalog.dat', usecols=(4,))          # velocities
dV = np.loadtxt('2768GCcatalog.dat', usecols=(5,))        # velocity errors

Rc = np.loadtxt('2768mag.dat', usecols=(4,))
Rcerr = np.loadtxt('2768mag.dat', usecols=(5,))
z = np.loadtxt('2768mag.dat', usecols=(6,))
zerr = np.loadtxt('2768mag.dat', usecols=(7,))

colour = Rc - z

n = len(RA)

#number of GC
print 'number of GCs:', n
GC = np.linspace(0,0,n)
ra = np.linspace(0,0,n)
dec = np.linspace(0,0,n)

for i in range(0,n):
    aux = SkyCoord(RA[i], DEC[i])             #reading RA and DEC from GC
    ra[i] = aux.ra.arcsec
    dec[i] = aux.dec.arcsec
    
        
y_rad = galcenter.dec.radian   # galaxy center DEC coordinate in radians
        
xmc = (ra - galcenter.ra.arcsec)
xm = xmc*(np.cos(y_rad))                       # setting the x coordinates of GC from its RA value relative to galactic center      
ym = (dec - galcenter.dec.arcsec)             # setting the y coordinates of the GC from its DEC value relative to galactic center
                     
             
xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # rotating coordinates 
ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    

cos_i = np.cos(ari_rad)               # cosine of inclination angle, used to shrink coordinates
sen_i = np.sin(ari_rad)

xsi = xs*cos_i
ysi = ys/cos_i

#colour bimodality

blue = []
red = []
blueGCx = []
blueGCy = []
redGCx = []
redGCy = []

blueGCx_plot = np.linspace(0,0,n)
blueGCy_plot = np.linspace(0,0,n)
redGCx_plot = np.linspace(0,0,n)
redGCy_plot = np.linspace(0,0,n)

c_sep = 0.57   #colour separation


for i in range(0,n):
    if colour[i] < c_sep:
        blue.append(V[i])
        blueGCx.append(xsi[i])
        blueGCy.append(ysi[i])
        blueGCx_plot[i] = (xs[i])
        blueGCy_plot[i] = (ys[i])
    else:
        red.append(V[i])
        redGCx.append(xsi[i])
        redGCy.append(ysi[i])
        redGCx_plot[i] = (xs[i])
        redGCy_plot[i] = (ys[i])     
        
blueGCx = np.asarray(blueGCx)
blueGCy = np.asarray(blueGCy)
redGCx = np.asarray(redGCx)
redGCy = np.asarray(redGCy)

print 'number of blue GC:', len(blue)
print 'number of red GC:', len(red)

 

r=np.sqrt(((xs**2)*cos_i)+((ys**2)/cos_i))   #distributing the GC along the radius
rblue = np.sqrt(((blueGCx**2)*cos_i)+((blueGCy**2)/cos_i))
rred = np.sqrt(((redGCx**2)*cos_i)+((redGCy**2)/cos_i))

cos_phi = xsi/r                              #azimuthal angle cosine & sine
sen_phi = ysi/r  

rvel = r

phi = np.arccos(cos_phi)               # azimuthal angle 

#checking for outliers

print r
for i in range(0, len(r)):
    if r[i] == np.amax(r):
        print 'color of outlier: ', colour[i]
        out_index = i

raux = np.delete(r, out_index)
diff = np.amax(r) - np.amax(raux)
if diff >= 100.00:
    r = raux

##############  binning section begins   ######################################################################



t = 0
n = r.size   #lenght of r --- number of objects to iterate on the following loops
for i in range(0, n):
    for j in range(i+1,n):
        if r[i] > r[j]:
               t = r[i]
               r[i] = r[j]
               r[j] = t
                  
nred = len(rred)
                  
for i in range(0, nred):
    for j in range(i+1,nred):
        if rred[i] > rred[j]:
               t = rred[i]
               rred[i] = rred[j]
               rred[j] = t
               
nblue = len(rblue)
               
for i in range(0, nblue):
    for j in range(i+1,nblue):
        if rblue[i] > rblue[j]:
               t = rblue[i]
               rblue[i] = rblue[j]
               rblue[j] = t
               
                      
                             
nbin = 6                     #number of bins
rgal = np.linspace(0,0,nbin+1)       #rgal contains the limits of each bin
rgalred = np.linspace(0,0,nbin+1)
rgalblue = np.linspace(0,0,nbin+1)


for h in range(1, nbin):
    rgal[h] = r[int(n*h/nbin)]     #binning
rgal[nbin] = r[n-1]
rgal[0] = r[0]


NGC = np.zeros(nbin)    #number of GC in each bin
rmed = [[] for i in range(0,nbin)]
for h in range(0, nbin):                              
     for i in range(0, n):
         if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
             NGC[h] = NGC[h] + 1
             rmed[h].append(r[i])
             
NGCreds = np.zeros(nbin)    #number of GC in each bin

for h in range(0, nbin):                              
     for i in range(0, nred):
         if (rgal[h+1] >= rred[i] and rred[i] > rgal[h]):
             NGCreds[h] = NGCreds[h] + 1

NGCblues = np.zeros(nbin)    #number of GC in each bin

for h in range(0, nbin):                              
     for i in range(0, nblue):
         if (rgal[h+1] >= rblue[i] and rblue[i] > rgal[h]):
             NGCblues[h] = NGCblues[h] + 1
             
#median value of the bins
radius = [[] for i in range(0,nbin)]
for i in range(0,nbin):
    radius[i] = np.linspace(rgal[i], rgal[i+1], 1000)   
median = np.linspace(0,0,nbin)
medianGC = np.linspace(0,0,nbin)
for i in range(0,nbin):
    median[i] = np.median(radius[i])
    medianGC[i] = np.median(rmed[i])

area = []
for h in range(0, nbin):
    area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
    
binsize = np.linspace(0,0,nbin)    
for h in range(0,nbin):
    binsize[h] = rgal[h+1]-rgal[h]  
binsize = binsize/2
    
dens2 = np.linspace(0,0,nbin)
densred = np.linspace(0,0,nbin)
densblue = np.linspace(0,0,nbin)
for h in range(0, nbin):
    dens2[h] = NGC[h]/area[h]
    densred[h] = NGCreds[h]/area[h]
    densblue[h] = NGCblues[h]/area[h]
    
#error in density

poi_err = (np.sqrt(NGC)/area)
dens2err = unp.uarray(dens2, poi_err)    
denslog2 = -2.5*np.log10(dens2)
denslog2err = -2.5*unp.log10(dens2err)  
poi_err2 = unp.std_devs(denslog2err)

poi_err_reds = (np.sqrt(NGCreds)/area)
densrederr = unp.uarray(densred, poi_err_reds)    
denslogred = -2.5*np.log10(densred)
denslogrederr = -2.5*unp.log10(densrederr)  
poi_errred = unp.std_devs(denslogrederr)

poi_err_blues = (np.sqrt(NGCblues)/area)
densblueerr = unp.uarray(densblue, poi_err_blues)    
denslogblue = -2.5*np.log10(densblue)
#for i in range(0,nbin):
#    if densblueerr[i] <= 0.0:
#        densblueerr[i] = 1
denslogblueerr = -2.5*unp.log10(densblueerr)
poi_errblue = unp.std_devs(denslogblueerr)


with open('GCdensity-circularbins.dat', 'w') as GC:        #write out the densities
    print >>GC, '# nbin   dens denslog err area NGC median rgal binsize'
    for h in range(0, nbin):
        print >>GC, h, dens2[h], denslog2[h], poi_err2[h], area[h], NGC[h], median[h], rgal[h+1], binsize[h]
with open('GCdensity-REDS.dat', 'w') as GC:        #write out the densities
    print >>GC, '# nbin   dens denslog err area NGC median rgal binsize'
    for h in range(0, nbin):
        print >>GC, h, densred[h], denslogred[h], poi_errred[h], area[h], NGCreds[h], median[h], rgal[h+1], binsize[h]
with open('GCdensity-BLUES.dat', 'w') as GC:        #write out the densities
    print >>GC, '# nbin   dens denslog err area NGC median rgal binsize' 
    for h in range(0, nbin):
        print >>GC, h, densblue[h], denslogblue[h], poi_errblue[h], area[h], NGCblues[h], median[h], rgal[h+1], binsize[h]


with open('GClog-circular.dat', 'w') as catalog:    
    for i in range(0,n+1):
        print >>catalog, phi[i]


##############################  plots #########################################
#some test plots to understand better what i'm doing in the end

#Colour distribution
hist = pylab.subplot(111)
hist.set_title('NGC2768 GC Color')
pylab.hist(colour, histtype='stepfilled', bins=20, color='gray')
pylab.plot((0.57, 0.57),(0.0, 25.0), 'k--')
hist.set_xlim([0.1, 1.2])
hist.set_ylim([0.0, 25.0])
pylab.xlabel('(Rc-z)')
pylab.ylabel('Number of GC')
pylab.show()
            
##color coding the GC with respect to colour (Rc - z)
#pylab.plot(-blueGCx_plot, blueGCy_plot, 'bo')
#pylab.plot(-redGCx_plot, redGCy_plot, 'ro')
#pylab.plot(0,0,'rx', markersize=8, label='galaxy center')
#pylab.xlabel('x (arcsec)')
#pylab.ylabel('y (arcsec)')                        
#pylab.legend(loc='lower left')
##pylab.show()          
#
##velocity distribution
#pylab.hist(V, bins=100)
#pylab.xlabel('Velocity')
#pylab.ylabel('Number of GC')
##pylab.show()
#
##pylab.errorbar(rmedian, dens, yerr=poi_err, marker='o', color='blue', linestyle='none', mfc='None' )
#pylab.plot(medianGC, denslog2, marker='D', color='green', linestyle='none', label='No Correction')
#pylab.xlabel('Radius')
#pylab.ylabel('GC Density')
#pylab.legend(loc='lower left')
#pylab.plt.gca().invert_yaxis()
##pylab.show()    
#        
#pylab.plot(rvel, V, 'ro')
#pylab.xlabel('R(arcseconds)')
#pylab.ylabel('V (km/s)')
#pylab.plt.gca().invert_yaxis()
##pylab.show()        

