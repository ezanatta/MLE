# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:42:22 2017
Code to plot metallicity vs probabilities and similar 

@author: emilio
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/home/emilio/MLE/2.2/')
from bin_functions import binning
from astropy import units as u
from astropy.coordinates import SkyCoord

def pne_radius(RA, DEC, gc, pa, i):
     
     pa = pa
     pa_rad = pa.to(u.rad)
     pa_rad = pa_rad.value  
     
     ra = np.linspace(0,0,len(RA))
     dec = np.linspace(0,0,len(RA))


     cos_i = np.cos(i.to(u.rad))
     yrad = np.cos(gc.dec.radian)  
     
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
     #ysi = ysi/np.sqrt(cos_i)
     
     return xsi, -ysi

gal = raw_input('Type Galaxy number:')
galinput = 'N'+gal+'input.dat'
galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput
galcat = '/home/emilio/MLE/2.5/N'+gal+'/prob_zh.dat'

RA = np.loadtxt(galcat, usecols=(0,))
DEC = np.loadtxt(galcat, usecols=(1,))
CaT = np.loadtxt(galcat, usecols=(2,))
CaTerr = np.loadtxt(galcat, usecols=(3,))
Z = np.loadtxt(galcat, usecols=(5,))
Zerr = np.loadtxt(galcat, usecols=(6,))
fb = np.loadtxt(galcat, usecols=(8,))
col = np.loadtxt(galcat, usecols=(9,))
like = np.loadtxt(galcat, usecols=(10,))
        
frej = []
zrej = []
colrej = []
        
for i in range(0, len(like)):
    if like[i] == 5.00:
         frej.append(fb[i])    
         zrej.append(Z[i])    
         colrej.append(col[i])

frej = np.asarray(frej)
zrej = np.asarray(zrej)
colrej = np.asarray(colrej)

plt.plot(fb, Z, marker='o', linestyle='None', color='black', alpha=0.8)
plt.plot(frej, zrej, marker='x',markersize=10, linestyle='None', color='red')
plt.title('NGC'+gal)
plt.xlabel('$L_{sph}(v, f)$', fontsize=20)
plt.ylabel('Z/H (sun)', fontsize=20)
plt.show()

plt.plot(Z, col, 'bo')
plt.xlabel('Z')
plt.ylabel('Color')
plt.show()
###################

galinput = 'N'+gal+'input.dat'
galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput

with open(galinput, 'r') as f:
    inp = [x.split(' ')[0] for x in open(galinput).readlines()]
        
c_sep = float(inp[4])
    
RAgal = inp[0]
DECgal = inp[1]
incl = float(inp[2])*u.deg
pa = float(inp[3])*u.deg
Vsys = float(inp[10])
d = float(inp[11])*u.mpc
galcenter = SkyCoord(RAgal, DECgal)

def atokpc(R, d):
    Dkpc = (d.value*np.pi*R)/648
    return Dkpc 

tab = np.genfromtxt('/home/emilio/MLE/2.5/N'+gal+'/prob_metal.dat', unpack=True)
rac = tab[0]
decc = tab[1]
xs = tab[2]
ys = tab[3]
v = tab[4]
fbt = tab[5]
liket = tab[7]
zt = tab[11]

xs, ys = pne_radius(rac, decc, galcenter, pa, incl) 

xs = atokpc(xs, d)
ys = atokpc(ys, d)

ztrej = []
xsrej = []
ysrej = []
vrej = []

for i in range(0, len(liket)):
    if liket[i] == 5.00:
        ztrej.append(zt[i])
        xsrej.append(xs[i])
        ysrej.append(ys[i])
        vrej.append(v[i])
        
ztrej = np.asarray(ztrej)        
xsrej = np.asarray(xsrej)
ysrej = np.asarray(ysrej)
vrej = np.asarray(vrej)

#plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
#        plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
#        plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
plt.plot(np.sqrt(xsrej**2+ysrej**2), vrej-Vsys, markersize=4, mec='red', mew=1, marker='x', mfc='None',linestyle='None', label='Rejected GC')
plt.scatter(np.sqrt(xs**2+ys**2), v-Vsys, c=zt, s=60, cmap='binary', alpha=0.8, vmin=min(zt)+1, vmax=max(zt)-1)
plt.xlim(0, 20)
plt.ylabel('$\Delta$V $(km/s)$', fontsize=20)
plt.xlabel('R $(Kpc)$', fontsize=20)
cb1 = plt.colorbar(fraction=0.04, orientation='vertical')
cb1.set_label('$Z/H (sun)$', fontsize=20)
#plt.plot((0, xmax),(0, 0), ':', label='Galaxy Systemic Velocity', color='black')
#plt.legend(loc='upper right', prop={'size':10})
plt.show() 

#plt.plot(np.sqrt(xs**2+ys**2), vfb, marker='o', color='purple', linestyle='None', label='GC')
#        plt.plot(rred, Vred, marker='o', color='red', linestyle='None')
#        plt.plot(rblue, Vblue, marker='o', color='blue', linestyle='None')
#plt.plot(np.sqrt(xsrej**2+ysrej**2), vrej-Vsys, markersize=4, mec='red', mew=1, marker='x', mfc='None',linestyle='None', label='Rejected GC')
plt.scatter(np.sqrt(xsrej**2+ysrej**2), vrej-Vsys, c=ztrej, s=60, cmap='binary', alpha=0.8, vmin=min(zt)+1, vmax=max(zt)-1)
plt.xlim(0, 20)
plt.ylabel('$\Delta$V $(km/s)$', fontsize=20)
plt.xlabel('R $(Kpc)$', fontsize=20)
cb1 = plt.colorbar(fraction=0.04, orientation='vertical')
cb1.set_label('$Z/H (sun)$', fontsize=20)
#plt.plot((0, xmax),(0, 0), ':', label='Galaxy Systemic Velocity', color='black')
plt.title('Only Rejected GCs - NGC '+gal)
plt.show() 

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(np.sqrt(xsrej**2+ysrej**2),vrej, ztrej, markersize=10, mec='red', mew=1, marker='x', mfc='None',linestyle='None', label='Rejected GC')
ax.scatter(np.sqrt(xs**2+ys**2),v, zt, s=50)
ax.set_xlim(0, 20)
ax.set_ylabel('$\Delta$V $(km/s)$', fontsize=20)
ax.set_xlabel('R $(Kpc)$', fontsize=20)
ax.set_zlabel('Z/H (sun)', fontsize=20)
plt.show() 

plt.plot(np.sqrt(xs**2+ys**2),zt, marker='o', linestyle='None', color='black', alpha=0.8)
plt.plot(np.sqrt(xsrej**2+ysrej**2), ztrej, marker='x',markersize=10, linestyle='None', color='red')
plt.title('NGC'+gal)
plt.xlabel('$R(kpc)$', fontsize=20)
plt.ylabel('Z/H (sun)', fontsize=20)
plt.show()

plt.scatter(zt, np.sqrt(xs**2+ys**2), c=fbt, s=50, cmap='Greys')
plt.plot(ztrej, np.sqrt(xsrej**2+ysrej**2), marker='x', markersize=10, linestyle='None', color='red')
plt.xlabel('Z/H (sun)', fontsize=20)
plt.ylabel('$R(Kpc)$', fontsize=20)
cb = plt.colorbar(fraction=0.04, orientation='vertical')
cb.set_label('$L_{sph}(v, f)$', fontsize=20)
plt.show()