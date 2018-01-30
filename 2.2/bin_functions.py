def binning(RA, DEC, V, galcenter, pa, i, nbins, op):
     import numpy as np
     from astropy import units as u
     from astropy.coordinates import SkyCoord
     
     pa = pa+90.00*u.deg
     pa_rad = pa.to(u.rad)
     pa_rad = pa_rad.value

     i_rad = i.to(u.rad)
     i_rad = i_rad.value

     n = len(RA)

     ra = np.linspace(0,0,n)
     dec = np.linspace(0,0,n)

     for i in range(0,n):
	 aux = SkyCoord(RA[i], DEC[i])             #reading RA and DEC from GC
	 ra[i] = aux.ra.arcsec
	 dec[i] = aux.dec.arcsec
       
   
   
     rad = np.zeros(n)  
     for i in range(0, n):
         gc = SkyCoord(RA[i], DEC[i])
         rad[i] = SkyCoord.separation(gc, galcenter).arcsec
     
     with open('checkradius-true.dat', 'w') as chkr:
         for k in range(0, n):
             print >>chkr, RA[k].value, DEC[k].value, rad[k]/60    
	
     y_rad = galcenter.dec.radian   # galaxy center DEC coordinate in radians
     #print 'Y RAD', y_rad
        
     xmc = (ra - galcenter.ra.arcsec)
     xm = xmc*(np.cos(y_rad))                       # trasnform the coordinates from spherical to planar   
     ym = (dec - galcenter.dec.arcsec)             # setting the y coordinates of the GC from its DEC value relative to galactic center
	
     xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # 	rotating coordinates 
     ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    
	
     cos_i = np.cos(i_rad)               # cosine of inclination angle, used to shrink coordinates
	
 
     xd = ((xs)*np.sqrt(cos_i))
     yd = ((ys)/np.sqrt(cos_i))
     
     r=np.sqrt((xs**2)*cos_i+((ys**2)/cos_i))   #distributing the GC along the radius
     
     
     
#     plt.plot(xs, ys, 'bo')
#     plt.plot(xm, ym, 'ro')
#     plt.plot(xs, yd, 'go')
#     plt.show()

     with open('checkradius.dat', 'w') as chkr2:
         for k in range(0, n):
             print >>chkr2, xd[k], yd[k], r[k]/60
	
     ###binning ####
     t = 0
     l = 0
     n = r.size   #lenght of r --- number of objects to iterate on the following loops
     for i in range(0, n):
         for j in range(i+1,n):
             if r[i] > r[j]:
               			t = r[i]
               			r[i] = r[j] #ordering r and V
               			r[j] = t
				l = V[i]
				V[i] = V[j]
				V[j] = l
	
     rgal = np.linspace(0,0,nbins+1)       #rgal contains the limits of each bin
     
     if op=='N':
         for h in range(1, nbins):
             rgal[h] = r[int(n*h/nbins)]     #binning
             rgal[nbins] = r[n-1]
             rgal[0] = r[0]     
     else:
         k = max(r)/nbins #size of each even spaced bins
         for h in range(1, nbins):
             rgal[h] = (h*k)
             rgal[nbins] = r[n-1]
             rgal[0] = r[0]  
	
     return r, rgal, V 
 
def density(nbin, r, rgal):
    import numpy as np
#    import uncertainties.unumpy as unp
    
    n = len(r)
    NGC = np.zeros(nbin)    #number of GC in each bin
    rmed = [[] for i in range(0,nbin)]
    for h in range(0, nbin):                              
         for i in range(0, n):
             if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
                 NGC[h] = NGC[h] + 1
                 rmed[h].append(r[i])
             
    #median value of the bins
    radius = [[] for i in range(0,nbin)]
    for i in range(0,nbin):
            radius[i] = np.linspace(rgal[i], rgal[i+1], 1000)   

    median = np.linspace(0,0,nbin)
    medianGC = np.linspace(0,0,nbin)
    for i in range(0,nbin):
        median[i] = np.median(radius[i])
        medianGC[i] = np.median(rmed[i])

    area = []
    for h in range(0, nbin):
        area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
    
    binsize = np.linspace(0,0,nbin)    
    for h in range(0,nbin):
        binsize[h] = rgal[h+1]-rgal[h]  
    binsize = binsize/2
    
    dens2 = np.linspace(0,0,nbin)
    for h in range(0, nbin):
        dens2[h] = NGC[h]/area[h]

    #error in density

    poi_err = np.sqrt(NGC)/area
    denslog2 = np.log10(dens2)    
#    dens2err = unp.uarray(dens2, poi_err)
#    denslog2err = unp.log10(dens2err)  
#    poi_err2 = unp.std_devs(denslog2err)
    
# manual error
    
#    poi_err2 = poi_err/(2.3*(dens2))
    poi_err2 = np.log10(poi_err)
    
    return dens2, denslog2, poi_err2, area, NGC, median, binsize, poi_err 
    
    
def bin_areafix(RA, DEC, V, galcenter, pa, i, nbins, op):         #Work in progress
     import numpy as np
     from astropy.coordinates import SkyCoord

     n = len(RA)

     r = np.zeros(n)     
     
     for i in range(0, n):
         gc = SkyCoord(RA[i], DEC[i])
         r[i] = SkyCoord.separation(gc, galcenter).arcsec
     
     with open('checkradius.dat', 'w') as chkr:
         for k in range(0, n):
             print >>chkr, RA[k].value, DEC[k].value, r[k]/60    
	
     ###binning ####
     t = 0
     l = 0
     n = r.size   #lenght of r --- number of objects to iterate on the following loops
     for i in range(0, n):
         for j in range(i+1,n):
             if r[i] > r[j]:
               			t = r[i]
               			r[i] = r[j] #ordering r and V
               			r[j] = t
				l = V[i]
				V[i] = V[j]
				V[j] = l
 
     rgal = np.linspace(0,0,nbins+1)       #rgal contains the limits of each bin
     
     ## check if there are isolated GCs
          
     n = len(r)
     if op=='N':
         for h in range(1, nbins):
             rgal[h] = r[int(n*h/nbins)]     #binning
             rgal[nbins] = r[n-1]
             rgal[0] = r[0]     
     else:
         k = max(r)/nbins #size of each even spaced bins
         for h in range(1, nbins):
             rgal[h] = (h*k)
             rgal[nbins] = r[n-1]
             rgal[0] = r[0]                            
                                 	
     return r, rgal, V   
     
def galaxy_inputs(gal):
    galinput = 'N'+gal+'input.dat'
    galinput = '../Galaxies/'+gal+'/'+galinput
    
    return galinput
    
def plotgcbins():

    from matplotlib import pyplot as plt
    from matplotlib.patches import Circle
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy import wcs
    from astropy.coordinates import SkyCoord
    import numpy as np
        
    def read_catalog_gc(gccat):
        RAg = np.loadtxt(gccat, usecols=(1,))
        DECg = np.loadtxt(gccat, usecols=(2,))
    
        return RAg, DECg

    with open('current_galaxy.dat', 'r') as f:     
        gal = [line.split()[0] for line in f]
        gal = gal[0]
    gc = '/home/emilio/MLE/Galaxies/'+gal+'/'+'N'+gal+'GC.dat'    
    galinput = 'N'+gal+'input.dat'
    galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput
    
    RA, DEC = read_catalog_gc(gc)
    coords = zip(RA, DEC)
    ra = []
    dec = [] 
    
    fits_file = '/home/emilio/MLE/Galaxies/'+gal+'/fits-null.fits'
    hdu = fits.open(fits_file)[0]
    hwcs = WCS(hdu.header)
    
    pix = hwcs.wcs_world2pix(coords, 1)

    for i in range(0,len(pix)):
        ra.append(pix[i][0])
        dec.append(pix[i][1])
        
    GCfile = '/home/emilio/MLE/2.2/test_GCdensity-circularbins.dat'

    r = np.loadtxt(GCfile, usecols=(7,))

    galcenter = SkyCoord.from_name('NGC '+gal)
    center = wcs.utils.skycoord_to_pixel(galcenter, hwcs)
    center = (float(center[0]), float(center[1]))        
        
    allgc = plt.subplot(1,1,1, projection=hwcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r') #you can change the contrast in the cmap options.
    plt.plot(ra, dec, marker='x',markerfacecolor='None', linestyle='none', markeredgewidth=1, markeredgecolor='blue', label='GCs')
    for r in r:        
        cir = Circle(center, r, fill=False)
        allgc.add_artist(cir)
    plt.legend(loc='lower right')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.show()
    
def plotgcbins_N(i):
    import numpy as np            
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    
    def plot_ellipse(semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e4,ax=None,plot_kwargs=None,fill=False,fill_kwargs=None):
        '''
                - create an ellipse in polar coordinates then transform to cartesian
                - if given an axes, plot an ellipse with plot_kwargs
                - if not given an axes, create a basic figure and axes and plot
                major keywords are:
                semimaj : length of semimajor axis
                semimin : length of semiminor axis
                phi : angle in radians semimajor axis is above positive x axis
                x_cent : x center
                y_cent : y center
                theta_num : number of points to sample from 0 to 2pi
        '''
        # Generate data for ellipse structure
        theta = np.linspace(0,2*np.pi,theta_num)
        r = semimaj*semimin / np.sqrt((semimin*np.cos(theta))**2 + (semimaj*np.sin(theta))**2)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        data = np.array([x,y])
        R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
        data = np.dot(R,data)
        data[0] += x_cent
        data[1] += y_cent

        # Plot!
        if ax == None:
                fig,ax = plt.subplots()

        if plot_kwargs == None:
                ax.plot(data[0],data[1],color='black',linestyle='-')        
        else:
                ax.plot(data[0],data[1],**plot_kwargs)


        if fill == True:
                ax.fill(data[0],data[1],**fill_kwargs)

        
    def read_chkr(gccat):
        x = np.loadtxt(gccat, usecols=(0,))
        y = np.loadtxt(gccat, usecols=(1,))
    
        return x, y
                 
    x, y = read_chkr('/home/emilio/MLE/2.2/checkradius.dat')
    r = np.loadtxt('/home/emilio/MLE/2.2/test_GCdensity-circularbins.dat', usecols=(7,))
    
    allgc = plt.subplot(1,1,1)
    plt.plot(x, y, marker='x',markerfacecolor='None', linestyle='none', markeredgewidth=1, markeredgecolor='blue', label='GCs')
    for r in r:        
        cir = Circle((0,0), r, fill=False)
        allgc.add_artist(cir)
    plt.legend(loc='lower right')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()    
    
def color_dist(colour):
    import numpy as np
    import matplotlib.pyplot as plt

    r = np.loadtxt('/home/emilio/MLE/2.2/test_GCdensity-circularbins.dat', usecols=(7,))   
    
    plt.plot(colour, r, 'bo')
    plt.show()

def bin_kpc(RA, DEC, RAgal, DECgal, d, nbins):
      from astropy.coordinates import SkyCoord
      
      n = len(RA)
      
      dist = list()
    
      g = SkyCoord(RAgal, DECgal, distance=d)
      
      for i in range(0, n):
          c = SkyCoord(RA[i], DEC[i], distance=d)
          sep = g.separation_3d(c)
          dist.append(sep)
      
      ###binning ####
      r = sorted(dist)
      n = len(r)   #lenght of r --- number of objects to iterate on the following loops
      
      rgal = list()      #rgal contains the limits of each bin
      for h in range(1, nbins):
        rgal.append(r[int(n*h/nbins)])     #binning
      rgal.insert(nbins, r[n-1])
      rgal.insert(0, r[0])
      
      return r, rgal 
    
def read_gc(gal):
    import numpy as np
    from astropy import units as u  
    
    galcat = 'N'+gal+'GC.dat' 
    galcat = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcat
    galinput = 'N'+gal+'input.dat'
    galinput = '/home/emilio/MLE/Galaxies/'+gal+'/'+galinput
    

    inp = [x.split(' ')[0] for x in open(galinput).readlines()]
    
    c_sep = float(inp[4])    
    RAgal = inp[0]
    DECgal = inp[1]
    i = float(inp[2])*u.deg
    pa = float(inp[3])*u.deg
    d = float(inp[11])*u.mpc
    #re = float(inp[12])
    
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    RA = RA*u.deg
    DEC = DEC*u.deg
    
    #inps = [RAgal, DECgal, i, pa, d, c_sep, re]
    inps = [RAgal, DECgal, i, pa, d, c_sep]
    return RA, DEC, inps
    
    
    
    
    
    
    
    
    
    
    