f = open('PNE', 'r')

with open('PNE', 'r') as f:
    RA = [line.split()[1] for line in f]
with open('PNE', 'r') as f:    
    DEC = [line.split()[2] for line in f]

print len(RA)
print len(DEC)

w = open('PNEcat.reg', 'w')

for h in range(0, 158):
    print >>w, RA[h], DEC[h]