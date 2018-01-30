# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 09:54:55 2017

@author: emilio
"""

#code only to plot N5866 GC number density - what a waste!

import numpy as np
import pylab
import uncertainties.unumpy as unp
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker

def plotgalfit(x, mue, mud, mud2, bn, re, n, h, h2):
    re_min = re/60
    hmin = h/60
    xmin = x/60
    srsclog = mue*np.exp(-bn*((xmin/re_min)**(1/n)-1))  
    exp1log = mud*np.exp(-xmin/hmin)
    exp2log = mud*np.exp(-xmin/hmin)
    totallog = srsclog+exp1log+exp2log   
 
    xlog = np.log10(xmin)

    return xlog, srsclog, exp1log, exp2log, totallog 
    
def galaxy_inputs(gal):
    galinput = 'N'+gal+'input.dat'
    galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput
    
    return galinput
    
def config_pne(N, area, r):
    area = area/3600
    
    denslog = N/area
    #denslog = np.log10(density)
    
    errlog = np.sqrt(N)/area
    #errlog = err/(2.3*(density))
    #errlog = np.log10(err)
    
    denslog10 = unp.uarray(denslog, errlog)
    denslog10 = unp.log10(denslog10)
    errlog = unp.std_devs(denslog10)
    
    rmin = r/60
        
    rlog = rmin
    #TO DO: add binsizes
    
    return denslog, errlog, rlog
    
def config_logplots(dens, err, r):
    
    denslog10 = unp.uarray(dens, err)
    denslog10 = unp.log10(denslog10)
    errlog = unp.std_devs(denslog10)
    
    #rmin = r/60
        
    rlog = r
    #TO DO: add binsizes

    return dens, errlog, rlog
    
    
galput = galaxy_inputs('5866')
inp = [x.split(' ')[0] for x in open(galput).readlines()]    
    
#n = 4.0
#bn = 2.0*n-0.327
#re = 12.9
#h = 30.19
#h2 = 13.65
#mue = 17.1
#mud = 16.1
#mud2 = 14.2

n = 4.0
bn = 2.0*n-0.327
re = 14.1
h = 34.29
h2 = 17.87
mue = 11.03
mud = 16.60
mud2 = 15.17

x = np.arange(0,1500,5.17)
nr = x.size

xlog, srsclog, exp1log, exp2log, totallog = plotgalfit(x, mue, mud, mud2, bn, re, n, h, h2)

rgc = np.loadtxt('density_from_rhodes.dat', usecols=(0,))
gcdens = np.loadtxt('density_from_rhodes.dat', usecols=(1,))
gcdens_err = np.loadtxt('density_from_rhodes.dat', usecols=(2,))

denslog, errlog, rlog = config_logplots(gcdens, gcdens_err, rgc)

#PNe
N = np.loadtxt('/home/emilio/MLE/2.2/PNEdensity-circularbins.dat', usecols=(5,))
A = np.loadtxt('/home/emilio/MLE/2.2/PNEdensity-circularbins.dat', usecols=(4,))
R = np.loadtxt('/home/emilio/MLE/2.2/PNEdensity-circularbins.dat', usecols=(7,))

denspne, errpne, rpne = config_pne(N, A, R)

plt = pylab.subplot(4,1,1)
plt.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
pylab.text(15,1.5, 'TOTAL')
pylab.plot((x/60), np.log10(totallog*3), color='black', linestyle='-', label = 'Total Light Profile (Spheroid+Disk)', dashes=[20,20,20,20,20,20])
pylab.errorbar(rlog, 2.5*np.log10(denslog), yerr=errlog, capsize=6, capthick=1, ecolor='magenta', marker = 'D', markeredgecolor='magenta', markersize=10, mfc='None', label='GC', linestyle='none')
pylab.errorbar(rpne, 2.5*np.log10(denspne), yerr=errpne, capsize=6, capthick=1, ecolor='green', marker = 'D', markeredgecolor='green', markersize=10, mfc='None', label='PNe', linestyle='none')
pylab.xlim(0.15575245742766164,17.911532604181087)
pylab.ylim(-5.7410379921259818,4.089095538057741)
         
plt2 = pylab.subplot(4,1,2)
plt2.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
pylab.text(15,1.5, 'SPHEROID')
pylab.plot((x/60), np.log10(srsclog*17), color='black', linestyle='--', label='Sersic Profile (Spheroid Light)')
pylab.errorbar(rlog, 2.5*np.log10(denslog), yerr=errlog, capsize=6, capthick=1, ecolor='magenta', marker = 'D', markeredgecolor='magenta', markersize=10, mfc='None', linestyle='none')
pylab.errorbar(rpne, 2.5*np.log10(denspne), yerr=errpne, capsize=6, capthick=1, ecolor='green', marker = 'D', markeredgecolor='green', markersize=10, mfc='None', label='PNe', linestyle='none')
pylab.xlim(0.15575245742766164,17.911532604181087)
pylab.ylim(-5.7410379921259818,4.089095538057741)

plt3 = pylab.subplot(4,1,3)
plt3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
plt3.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
plt3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
pylab.text(15,1.5, 'DISK 1')
pylab.plot((x/60), np.log10(exp1log*6), color ='black', linestyle='-', label='Exp Disk Profile (Disk Light)', dashes=[8, 4, 2, 4, 2, 4])                
pylab.errorbar(rlog, 2.5*np.log10(denslog), yerr=errlog, capsize=6, capthick=1, ecolor='magenta', marker = 'D', markeredgecolor='magenta', markersize=10, mfc='None', linestyle='none')
pylab.errorbar(rpne, 2.5*np.log10(denspne), yerr=errpne, capsize=6, capthick=1, ecolor='green', marker = 'D', markeredgecolor='green', markersize=10, mfc='None', label='PNe', linestyle='none')
pylab.xlim(0.15575245742766164,17.911532604181087)
pylab.ylim(-5.7410379921259818,4.089095538057741)

plt4 = pylab.subplot(4,1,4)
plt4.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
plt4.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
plt4.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
pylab.text(15,1.5, 'DISK 2')
pylab.plot((x/60), np.log10(exp2log*6), color ='black', linestyle='-', label='Exp Disk Profile (Disk Light)', dashes=[8, 4, 2, 4, 2, 4])                
pylab.errorbar(rlog, 2.5*np.log10(denslog), yerr=errlog, capsize=6, capthick=1, ecolor='magenta', marker = 'D', markeredgecolor='magenta', markersize=10, mfc='None', linestyle='none')
pylab.errorbar(rpne, 2.5*np.log10(denspne), yerr=errpne, capsize=6, capthick=1, ecolor='green', marker = 'D', markeredgecolor='green', markersize=10, mfc='None', label='PNe', linestyle='none')
pylab.xlim(0.15575245742766164,17.911532604181087)
pylab.ylim(-5.7410379921259818,4.089095538057741)
              
                                
pylab.subplots_adjust(hspace=.0001)
pylab.setp(plt.get_xticklabels(), visible=False)
pylab.setp(plt2.get_xticklabels(), visible=False)
plt.set_ylabel('log($\Sigma_{T}$)', fontsize=20)
plt2.set_ylabel('log($\Sigma_{S}$)', fontsize=20)
plt3.set_ylabel('log($\Sigma_{D}$)', fontsize=20)
plt4.set_ylabel('log($\Sigma_{D}$)', fontsize=20)
plt4.set_xlabel('R ($arcmin$)', fontsize=20)
pylab.show()

########## comparing power law and r1/4 fits

ap, ap1 = 0.96, -1.76
av, av1 = 3.52, -2.52

pvauc = np.exp(av+av1*(rlog**(1/4)))
plaw = np.exp(ap+ap1*np.log10(rlog))

plt5 = pylab.subplot(1, 1, 1)
plt5.plot(rlog, pvauc, linestyle='--', color='black')
plt5.plot(rlog, plaw, linestyle='-', color='black')
plt5.plot(rlog, denslog, 'bo')

 