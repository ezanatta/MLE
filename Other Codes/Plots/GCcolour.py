import numpy as np
import pylab

Rc = np.loadtxt('2768mag.dat', usecols=(4,))
print Rc
Rcerr = np.loadtxt('2768mag.dat', usecols=(5,))
z = np.loadtxt('2768mag.dat', usecols=(6,))
print z
zerr = np.loadtxt('2768mag.dat', usecols=(7,))

colour = Rc - z

#Colour distribution
pylab.hist(colour, bins=20)
pylab.xlabel('(Rc-z)')
pylab.ylabel('Number of GC')
pylab.show()