# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 17:32:21 2018

@author: emilio
"""

from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option('-f','--file', dest='filename', help='file name')
parser.add_option('-c','--column', dest='column', help='column to read fb')

options, args = parser.parse_args()


fb = np.loadtxt(args[0], usecols=(int(args[1]),))   

ndisk = 0
nsph = 0

for item in fb:
    if item >= 0.5:
        nsph = nsph+1
    else:
        ndisk = ndisk+1
        
print('Disk GCs = %d'%ndisk)
print('Sph GCs = %d'%nsph)
