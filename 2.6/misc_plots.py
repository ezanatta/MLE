# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:34:57 2017

Miscelaneous plots for Chromodynamical Analysis of SOs... (2017-2018)

@author: emilio
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
Throughout the code, the order of the galaxies is the following, in all arrays:

1023, 2768, 3115, 7457

"""


def ngc(galaxy):
    
    p_a_sph = sum(np.loadtxt('home/emilio/MLE/to_send_ari/prob_out'+galaxy+'.dat', usecols=(6,))
    p_a_disk = 1-p_r_sph
    clr = np.loadtxt('home/emilio/MLE/to_send_ari/prob_out'+galaxy+'.dat', usecols=(7,))
    
    c_sep = np.loadtxt('home/emilio/MLE/Galaxies/'+galaxy+'/N'+galaxy'+input.dat', usecols=(0,))[4]

    for i in range(0, len(clr)):
        if clr[i] >= c_sep:
            p_r_disk.append(p_a_disk[i])
        else:
            p_r_sph.append(p_a_sph[i])

    
            


            



#from Cortesi2013b, the Bulge/Total ratio and, from Alabi2017, the total mass for each galaxy (in log $M\odot$):

gal = ['N1023', 'N2768', 'N3115', 'N7457']

bt = np.array([0.53, 0.71, 0.74, 0.30])
mass = np.array([10.99,11.21,10.93,10.13])
mass_nolog = 10**mass
vs = np.array([5.3,3.3,3.3,3.7])
vserr_max = np.array([0.8,0.6,0.6,0.6])
vserr_min = np.array([0.7,0.6,0.7,0.7])

#number of GCs in the bulge and disk, separated by colour, for each galaxy:

dummy = np.array([np.nan, np.nan, np.nan, np.nan])


red_bulge = np.array([13,45,41,np.nan])
red_disk = np.array([19,12,40,np.nan])
blue_bulge = np.array([18.6,42,32,np.nan])
blue_disk = np.array([25.4,6,36,np.nan])

all_disk = np.array([44.4,18,76,34])
all_bulge = np.array([31.6,87,73,6])
rej = np.array([15,10,32,4])

nb = all_bulge/mass_nolog
nd = all_disk/mass_nolog
rb = red_bulge/mass_nolog
bb = blue_bulge/mass_nolog

mark = ['s', 'o', '^', 'D']
fill=['none', 'full', 'none', 'full']

fig, ax = plt.subplots()
ax.set_xlabel('B/T', fontsize=20)
ax.set_ylabel('$N_{GC}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(bt, nb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle='full', color='black', label=g, linestyle='none', markersize=10)
for x, y, m, f in zip(bt, nd, mark, fill):
    ax.plot([x], [y], marker=m, fillstyle='none', color='black', linestyle='none', markersize=10)
leg = ax.legend(numpoints=1, frameon=False, fontsize=10)
leg.legendHandles[0]._sizes = [20]
leg.legendHandles[1]._sizes = [20]
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('B/T', fontsize=20)
ax.set_ylabel('$N_{GC}^{Red, Sph}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(bt, rb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='red', label=g, linestyle='none', markersize=10)
leg = ax.legend(numpoints=1, frameon=False, fontsize=10, loc='upper left')
leg.legendHandles[0]._sizes = [20]
leg.legendHandles[1]._sizes = [20]
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('B/T', fontsize=20)
ax.set_ylabel('$N_{GC}^{Sph}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(bt, bb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='black', label=g, linestyle='none', markersize=10)
for x, y, m, f, g in zip(bt, bb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='blue', linestyle='none', markersize=10)
for x, y, m, f, g in zip(bt, rb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='red', linestyle='none', markersize=10)
leg = ax.legend(numpoints=1, frameon=False, fontsize=10, loc='upper left')
leg.legendHandles[0]._sizes = [20]
leg.legendHandles[1]._sizes = [20]
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('$V/\sigma$', fontsize=20)
ax.set_ylabel('$N_{GC}^{Red,Sph}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(vs, rb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='red', label=g, linestyle='none', markersize=10)   
#for x, y, m, f in zip(vs, nd, mark, fill):
#    ax.plot([x], [y], marker=m, fillstyle=f, color='red', linestyle='none')
leg = ax.legend(numpoints=1, frameon=False, fontsize=10)
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('$V/\sigma$', fontsize=20)
ax.set_ylabel('$N_{GC}^{Sph}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(vs, nb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='magenta', label=g, linestyle='none', markersize=10)   
#for x, y, m, f in zip(vs, nd, mark, fill):
#    ax.plot([x], [y], marker=m, fillstyle=f, color='red', linestyle='none')
leg = ax.legend(numpoints=1, frameon=False, fontsize=10)
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel('$V/\sigma$', fontsize=20)
ax.set_ylabel('$N_{GC}^{Disk}/M_{*}$', fontsize=20)
for x, y, m, f, g in zip(vs, nd, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle=f, color='green', label=g, linestyle='none', markersize=10)   
#for x, y, m, f in zip(vs, nd, mark, fill):
#    ax.plot([x], [y], marker=m, fillstyle=f, color='red', linestyle='none')
leg = ax.legend(numpoints=1, frameon=False, fontsize=10)
plt.show()

#WISE mag values for all galaxies (WISE public database)

w1 = np.array([7.822, 9.022, 6.625, 10.078])
w2 = np.array([7.870, 9.032, 7.350, 10.089])
w3 = np.array([7.346, 8.033, 6.664, 9.478])

rm = red_bulge/mass_nolog
rm[3] = 36/mass_nolog[3]
red_bulge[3] = 36

coly = w1-w2
colx = w2-w3

fig, ax = plt.subplots()
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='5%', pad=0.05)
ax.set_xlabel('[4.6]-[12] (mag)', fontsize=20)
ax.set_ylabel('[3.4]-[4.6] (mag)', fontsize=20)
ax.set_ylim(-1,4)
ax.set_xlim(-1,7)
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.xaxis.grid(color='gray', linestyle='dashed')
for x, y, m, f, g in zip(dummy, bb, mark, fill, gal):
    ax.plot([x], [y], marker=m, fillstyle='full', color='black', label=g, linestyle='none', markersize=10)
for x, y, r, m, f, g in zip(colx, coly, red_bulge, mark, fill, gal):
    ax.scatter([x], [y], c=r, marker=m, s=100, vmin=10, vmax=40)
#cb1 = fig.colorbar(, cax=cax, orientation='vertical')
#cb1.set_label('$N_{GC}^{red, sph}$', fontsize=20)
#for x, y, m, f in zip(vs, nd, mark, fill):
#    ax.plot([x], [y], marker=m, fillstyle=f, color='red', linestyle='none')
ax.legend(numpoints=1, frameon=False, fontsize=10)
plt.show()