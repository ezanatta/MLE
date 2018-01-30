# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:53:37 2017

@author: emilio
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

font = {'family': ['Computer Modern'],
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

def bi(op):
    
    if op == 'allGC':
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R   
        
    if op=='redGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
         
        
    if op=='blueGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
                
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax

def uni(op):
        
    V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
    Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
    Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
    sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
    sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
    sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
    R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
    Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
    Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
    Rmin = R - Rmin
    Rmax = Rmax - R
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax

def vsigma():   
    
        V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax = bi('allGC')
        V_r, Verr_r_min, Verr_r_max, sigma_r, sigmaerr_r_min, sigmaerr_r_max, R_r, Rmin_r, Rmax_r = bi('redGC')
        V_b, Verr_b_min, Verr_b_max, sigma_b, sigmaerr_b_min, sigmaerr_b_max, R_b, Rmin_b, Rmax_b = bi('blueGC')
   
        Vos = V/sigma
        Vosreds = V_r/sigma_r
        Vosblues = V_b/sigma_b
        
        Vos_max = Verr_max/sigmaerr_max
        Vos_min = Verr_min/sigmaerr_min
        Vosreds_max = Verr_r_max/sigmaerr_r_max
        Vosreds_min = Verr_r_min/sigmaerr_r_min
        Vosblues_max = Verr_b_max/sigmaerr_b_max
        Vosblues_min = Verr_b_min/sigmaerr_b_min
        
        return Vos, Vosreds, Vosblues, Vos_max, Vos_min, Vosreds_max, Vosreds_min, Vosblues_max, Vosblues_min, R, R_r, R_b
        
def vsigma_uni():   
    
        V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax = bi('allGC')
   
        Vos = V/sigma
        
        Vos_max = Verr_max/sigmaerr_max
        Vos_min = Verr_min/sigmaerr_min
        
        return Vos, Vos_max, Vos_min, R
        
def atokpc(R):
    
    galinput = 'N'+gal+'input.dat'
    galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput
    inp = [x.split(' ')[0] for x in open(galinput).readlines()]
    d = float(inp[11])
    Dkpc = (d*np.pi*R)/648
    return Dkpc   
    
#how to proceed with error propagation with assymetrical errorbars? let's proceed without it for now:
    
gals = ['2768', '3115']

Vs = []
Vsmax = []
Vsmin = []
Vsreds = []
Vsblues = []
Vsredsmax =  []
Vsredsmin = []
Vsbluesmax = []
Vsbluesmin = []
Rs = []
Rsreds = []
Rsblues = []

for gal in gals:
    Vos, Vosreds, Vosblues, Vos_max, Vos_min, Vosreds_max, Vosreds_min, Vosblues_max, Vosblues_min, R, R_r, R_b = vsigma()
    R = atokpc(R)
    R_r = atokpc(R_r)
    R_b = atokpc(R_b)
    
    Vs.append(Vos)
    Vsmax.append(Vos_max)
    Vsmin.append(Vos_min)
    Vsreds.append(Vosreds)
    Vsredsmax.append(Vosreds_max)
    Vsredsmin.append(Vosreds_min)
    Vsblues.append(Vosblues)
    Vsbluesmax.append(Vosblues_max)
    Vsbluesmin.append(Vosblues_min)
    Rs.append(R)
    Rsreds.append(R_r)
    Rsblues.append(R_b)
    
gal = '7457'

Vos, Vos_max, Vos_min, R = vsigma_uni()
R = atokpc(R)
    
Vs.append(Vos)
Vsmax.append(Vos_max)
Vsmin.append(Vos_min)
Rs.append(R)


f, ax = plt.subplots(3, 3, sharex=True,sharey=True)
         
for i in range(0, 3):   
    
    ax[i][0].errorbar(Rs[i], Vs[i], yerr=[Vsmin[i], Vsmax[i]], markerfacecolor='none', markersize=10, ecolor='purple', markeredgecolor='purple', marker='o', linestyle='none')
    ax[i][0].set_ylim(0.0, 5.0)

for i in range(0, 2):

    ax[i][1].errorbar(Rsreds[i], Vsreds[i], yerr=[Vsredsmin[i], Vsredsmax[i]], markerfacecolor='none',ecolor='red', markersize=10,markeredgecolor='red', marker='o', linestyle='none')
    
for i in range(0, 2):
    
    ax[i][2].errorbar(Rsblues[i], Vsblues[i], yerr=[Vsbluesmin[i], Vsbluesmax[i]], markerfacecolor='none',ecolor='blue', markersize=10,markeredgecolor='blue', marker='o', linestyle='none')
    
plt.setp(ax[0][0].get_xticklabels(), visible=False)  
plt.setp(ax[0][1].get_xticklabels(), visible=False)
plt.setp(ax[0][1].get_yticklabels(), visible=False) 
plt.setp(ax[1][1].get_yticklabels(), visible=False)
plt.setp(ax[1][2].get_yticklabels(), visible=False) 
plt.setp(ax[1][1].get_xticklabels(), visible=True)
plt.setp(ax[1][2].get_xticklabels(), visible=True)  
plt.setp(ax[0][2].get_xticklabels(), visible=False)
plt.setp(ax[0][2].get_yticklabels(), visible=False) 
loc = plticker.MultipleLocator(base = 1.1)
loc2 = plticker.MultipleLocator(base = 5.0)


ax[1][0].set_ylabel('$V/\sigma$', fontsize=20)

ax[1][2].set_xlabel('$R(kpc)$', fontsize=20)
ax[2][0].set_xlabel('$R(kpc)$', fontsize=20)
ax[1][1].set_xlabel('$R(kpc)$', fontsize=20)


for i in range(0, 3):
    #ax[0][i].set_ylim(-1, 4)
    ax[0][i].yaxis.set_major_locator(loc)
    ax[0][i].xaxis.set_major_locator(loc2)
    ax[0][i].set_xticks((5, 10, 15, 20))
    ax[0][i].minorticks_on()

plt.subplots_adjust(hspace=.0001)
plt.subplots_adjust(wspace=.0001)
f.delaxes(ax[2][1])
f.delaxes(ax[2][2])
plt.savefig('vsigma.png', dpi=300, format='png')
plt.show()
