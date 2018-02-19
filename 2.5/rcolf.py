# -*- coding: utf-8 -*-
"""
Created on Mon May 29 17:09:34 2017

@author: emilio
"""
def plots_from_prob():
    import matplotlib.pyplot as plt
    import numpy as np    
    
    f = np.loadtxt(f_file, usecols=(0,))
    fb = np.loadtxt('unliky.dat', usecols=(6,))
    xs = np.loadtxt('unliky.dat', usecols=(2,))
    ys = np.loadtxt('unliky.dat', usecols=(3,))
    col = np.loadtxt('unliky.dat', usecols=(9,))
    like = np.loadtxt('unliky.dat', usecols=(10,))
    
    xrej = []
    yrej = []
    colrej = []
    
    for i in range(0, len(like)):
        if np.abs(like[i]) > 0.00:
            xrej.append(xs[i])      
            yrej.append(ys[i])    
            colrej.append(col[i])            
    print xrej         
    
    fig = plt.subplot(1, 1, 1)
    plt.hist(f, bins=20, color='gray')
    fig.set_xlabel('$f_{i}$')
    fig.set_ylabel('N')
    plt.show()
    fig2 = plt.subplot(1, 1, 1)
    plt.hist(fb, bins=20, color='gray')
    fig2.set_xlabel('$L_{Sph}(v_{i}, f_{i})$')
    fig2.set_ylabel('N')
    plt.show()
    fig3 = plt.subplot(1, 1, 1)
    plt.scatter(col, np.sqrt(xs**2+ys**2), c=fb, cmap='Greys_r', vmin=0.0, vmax=1.0, s=50)
    for i in range(0, len(colrej)):
        fig3.plot(colrej[i], np.sqrt(xrej[i]**2+yrej[i]**2), 'ro',linestyle='None', label='Rejected GC')    
    fig3.set_xlabel('(B-R)(mag)')
    fig3.set_ylabel('R$(arcsec)$')
    cb = plt.colorbar(fraction=0.05)
    cb.set_label('$L_{Sph}(v_{i}, f_{i})$')
    plt.show()
    
#    plt.hist([f, fb], bins=20, histtype='barstacked', color=['gray', 'white'])
#    plt.xlabel('$f_{i}$')
#    plt.ylabel('N')    
#    plt.show()
    
def atokpc(R, d):
    Dkpc = (d.value*3.141592653589793*R)/648
    return Dkpc 
