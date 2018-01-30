def HMS2deg(ra='', dec=''):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    D, M, S = [float(i) for i in dec.split()]
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    H, M, S = [float(i) for i in ra.split()]
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC
    
import numpy as np    
    
print 'Input input catalog name: '
catname = raw_input()    
  
ras = []
decl = []    
    
with open(catname, 'r') as cat:
    ra = np.loadtxt('PNEcatalogwithspaces.dat', usecols=(0,))
    dec = np.loadtxt('PNEcatalogwithspaces.dat', usecols=(1,))        
    for i in range(0,len(ra)):
        ras.append(HMS2deg(ra=ra[i]))
        decl.append(HMS2deg(dec=dec[i]))
        print i
        
print 'Input output catalog name: '
out = raw_input()
    
with open(out, 'w') as o:
    for i in range(0,len(ra)):
        print >>o, ras[i], ' ', decl[i]
