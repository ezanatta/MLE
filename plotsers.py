# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:35:59 2017

@author: emilio
"""

import numpy as np
import matplotlib.pyplot as plt

def sersic(x, n):
    bn = 2*n-0.237
    s = 20 + (2.5*bn/np.log(10))*((x**1/n)-1)
    #s = 10e-8*np.exp(-bn*((x**(1/n))-1))    
    return s

n = [0.5, 1, 4, 7, 10]    
c = ['blue', 'red', 'black', 'magenta', 'green' ]   

x = np.arange(0, 10, 0.3) 

for i in range(0, len(n)):  
    plt.plot(x, sersic(x, n[i]), color=c[i], label='n = '+str(n[i]))
plt.gca().invert_yaxis() 
plt.legend()    
plt.show()