import numpy as np

f = open('PNE', 'r')

with open('PNevel', 'r') as f:
    RA = [line.split()[1] for line in f]
with open('PNevel', 'r') as f:    
    DEC = [line.split()[2] for line in f]



RAlist = np.linspace(0,0,len(RA))
DEClist = np.linspace(0,0,len(RA))

for i in range(0,len(RA)):
    RAaux = list(RA[i])
    DECaux = list(DEC[i])
    RA[i] = RAaux[0]+RAaux[1]+'h'+RAaux[3]+RAaux[4]+'m'+RAaux[6]+RAaux[7]+RAaux[8]+RAaux[9]+RAaux[10]+'s'
    DEC[i] = '+' + DECaux[0] + DECaux[1] + 'd' + DECaux[3] + DECaux[4] + 'm'+ DECaux[6] + DECaux[7] + DECaux[8]+ DECaux[9] +'s'
    
    
w = open('PNEcoords.dat', 'w')

for h in range(0, 158):
    print >>w, RA[h], DEC[h]