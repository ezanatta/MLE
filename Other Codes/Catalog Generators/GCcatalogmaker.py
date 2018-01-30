def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1}d{2}m{3}s'.format(ds, deg, decM, decS)
  
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1}h{2}m{3}s'.format(rs, raH, raM, raS)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC
    
import numpy as np    
    
print 'Input input catalog name: '
catname = raw_input()    
    
n = 108    
ras = []
decl = []    
    
with open(catname, 'r') as cat:
    ra = np.loadtxt('2768GCcatalog.dat', usecols=(2,))
    dec = np.loadtxt('2768GCcatalog.dat', usecols=(3,))        
    for i in range(0,len(ra)):
        ras.append(deg2HMS(ra=ra[i],round=True))
        decl.append(deg2HMS(dec=dec[i],round=True))
        print i
        
print 'Input output catalog name: '
out = raw_input()
    
with open(out, 'w') as o:
    for i in range(0,len(ra)):
        print >>o, ras[i], ' ', decl[i]