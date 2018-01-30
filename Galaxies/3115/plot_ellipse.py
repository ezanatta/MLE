import numpy as np
import pylab

#photometric parameters (fitted by galfit)
n = 4.0      #sersic index
bn = 2.0*n-0.327
re = 12.8952    #effective radius (pix)
h = 30.1882   #disk scale lenght (exp1)  (pix)
h2 = 13.6473    #disk scale lenght (exp2)  (pix)
mue = 17.08302541201160   
mud = 16.08248348110908
mud2 = 14.26003386396783

x = np.arange(0,800,5.17)
rad = np.loadtxt('model.dat', usecols=(1,))
rr = np.loadtxt('model.dat', usecols=(2,))
err = np.loadtxt('model.dat', usecols=(3,))
logrr = -2.5*np.log10(rr)+1
logerr = -2.5*0.4343*err/rr

dens = np.loadtxt('PNEdensity-circularbins.dat', usecols=(2,))
denserr = np.loadtxt('PNEdensity-circularbins.dat', usecols=(3,))
rmed = np.loadtxt('PNEdensity-circularbins.dat', usecols=(6,))
binsize = np.loadtxt('PNEdensity-circularbins.dat', usecols=(8,))

dens = dens-5.9
nr = x.size

total = np.linspace(0,0,nr)
srsc = np.linspace(0,0,nr)
exp1 = np.linspace(0,0,nr)
exp2 = np.linspace(0,0,nr)

for i in range(0,nr):
    srsc[i] = mue*np.exp(-bn*((x[i]/re)**(1/n)-1))
    exp1[i] = mud*np.exp(-x[i]/h)
    exp2[i] = mud2*np.exp(-x[i]/h2)
    total[i] = srsc[i]+exp1[i]+exp2[i]
    
total = 17.08*np.exp(-7.67*(((x/(12.9))**(1/4.))-1))+16.07*np.exp(-x/(30.19))+14.26*np.exp(-x/(13.65))   
total = -2.5*np.log10(total)
srsc = -2.5*np.log10(srsc)
exp1 = -2.5*np.log10(exp1)
exp2 = -2.5*np.log10(exp2)

pylab.plot(rad, logrr, marker='D', linestyle='none', markersize=5, color='red', label='Ellipse Fit', mfc='None')
pylab.plot(x, srsc, color='green', label='Sersic', linestyle= '--')
pylab.plot(x, exp1, color ='cyan', label='Exp Disk 1', linestyle= '--')
pylab.plot(x, exp2, color='red', label= 'Exp Disk 2', linestyle= '--')
pylab.plot(x, total, color='black', label = 'Total Profile', linewidth=1)
pylab.errorbar(rmed, dens, yerr=denserr, xerr=binsize, marker = 'o',color='green', markersize=8, label='PNe Density', linestyle='none', mfc='None')
pylab.xlabel('$R (arcsec)$', fontsize=14)
pylab.ylabel('$\Sigma$ ($N/arcsec^2$)', fontsize=14)
pylab.legend(loc='upper right')
pylab.xlim(0,650)
pylab.ylim(-5,15)
pylab.plt.gca().invert_yaxis()
pylab.show()