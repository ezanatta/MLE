# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:57:02 2017

@author: emilio
"""


import matplotlib.pyplot as plt
import numpy as np    

gal = '3115'    
    
f_file = '/home/emilio/MLE/Galaxies/'+gal+'/N'+gal+'_f.GC.dat'    
    
f = np.loadtxt(f_file, usecols=(0,))
fb = np.loadtxt('for_lod_0.dat', usecols=(6,))
xs = np.loadtxt('for_lod_0.dat', usecols=(2,))
ys = np.loadtxt('for_lod_0.dat', usecols=(3,))
col = np.loadtxt('for_lod_0.dat', usecols=(7,))
like = np.loadtxt('for_lod_0.dat', usecols=(8,))
v = np.loadtxt('for_lod_0.dat', usecols=(4,))
    
xrej = []
yrej = []
colrej = []
vrej = []
    
for i in range(0, len(like)):
    if like[i] == 5.00:
        xrej.append(xs[i])    
        yrej.append(ys[i])    
        colrej.append(col[i])
        vrej.append(v[i])            
print xrej         
    
fig = plt.subplot(1, 1, 1)
plt.hist(f, bins=20, color='gray')
fig.set_xlabel('$f_{i}$', fontsize=20)
fig.set_ylabel('N', fontsize=20)
plt.show()
    
fig2 = plt.subplot(1, 1, 1)
plt.hist(fb, bins=20, color='gray')
fig2.set_xlabel('$L_{Sph}(v_{i}, f_{i})$', fontsize=20)
fig2.set_ylabel('N', fontsize=20)
plt.show()
    
fig3 = plt.subplot(1, 1, 1)
plt.scatter(col, np.sqrt(xs**2+ys**2), c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0, s=50)
for i in range(0, len(colrej)):
    fig3.plot(colrej[i], np.sqrt(xrej[i]**2+yrej[i]**2), 'ro',linestyle='None', label='Rejected GC')    
fig3.set_xlabel('(g-i)(mag)', fontsize=20)
fig3.set_ylabel('R$(arcsec)$', fontsize=20)
cb = plt.colorbar(fraction=0.05)
cb.set_label('$L_{Sph}(v_{i}, f_{i})$', fontsize=20)
plt.show()

xrej = np.asarray(xrej)
yrej = np.asarray(yrej)
vrej = np.asarray(vrej)

plt.plot(np.sqrt(xs**2+ys**2), v, color='purple', linestyle='None')
plt.plot(np.sqrt(xrej**2+yrej**2), vrej, color='red', linestyle='None')
plt.show()

plt.plot(xs, v, marker='o',color='purple', linestyle='None')
plt.plot(xrej, vrej, marker='o',color='red', linestyle='None')
plt.show()
    
#    plt.hist([f, fb], bins=20, histtype='barstacked', color=['gray', 'white'])
#    plt.xlabel('$f_{i}$')
#    plt.ylabel('N')    
#    plt.show()
    