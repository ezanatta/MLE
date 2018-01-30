from astropy.coordinates import SkyCoord
import numpy as np
import os

print('This is a catalog maker for PNs catalogs fit into MLE codes.')

gal = raw_input('Enter galaxy number: ')
cat_input = raw_input('Enter the input catalog: ')



cat_input = '/home/emilio/MLE/Galaxies/'+gal+'/'+cat_input

lines = []

rah = []
ram = []
ras = []
decd = []
decm = []
decs = []

ra = []
dec = []

V = np.loadtxt(cat_input, usecols=(6,))

f = open(cat_input, 'r')

for column in f:
    lines.append(column)
    
for i in range(0,len(lines)-1):
    print i
    rah.append(lines[i][1]+lines[i][2]+'h')
    ram.append(lines[i][4]+lines[i][5]+'m')
    ras.append(lines[i][7]+lines[i][8]+lines[i][9]+lines[i][10]+'s')
    decd.append(lines[i][12]+lines[i][14]+lines[i][15]+'d')
    decm.append(lines[i][17]+lines[i][18]+'m')
    decs.append(lines[i][20]+lines[i][21]+lines[i][22]+lines[i][23]+'s')
    
for i in range(0,len(lines)-1):
    ra.append(rah[i]+ram[i]+ras[i])
    dec.append(decd[i]+decm[i]+decs[i])

    
with open('intermediate_catalog.dat', 'w') as o:
    for i in range(0,len(ra)):
        print >>o, ra[i], dec[i]

with open('intermediate_catalog.dat', 'r') as f:
    RA = [line.split()[0] for line in f]              #PNe catalog
with open('intermediate_catalog.dat', 'r') as f:    
    DEC = [line.split()[1] for line in f]
    
for i in range(0, len(RA)):
    c = SkyCoord(RA[i], DEC[i])
    RA[i], DEC[i] = c.ra.degree, c.dec.degree

cat_out = 'N'+gal+'PNE.dat'
cat_out = '/home/emilio/MLE/Galaxies/'+gal+'/'+cat_out
    
with open(cat_out,'w') as pnecat:
    for i in range(0, len(RA)):
        print >>pnecat, i+1, "%.8f"%RA[i], "%.8f"%DEC[i], V[i]

f.close()
os.system('rm intermediate_catalog.dat')
    
