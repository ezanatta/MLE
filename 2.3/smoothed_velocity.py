import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
import sys
sys.path.insert(0, '../2.2/')
from bin_functions import galaxy_inputs
from scipy import stats
from astropy.io import fits
from astropy.wcs import WCS

def plot_vmap_singleplot(path, vn):
    
    hdu_list = fits.open(path+'fits-null.fits')
    hdu = hdu_list[0]
    wcs = WCS(hdu.header)
    
    DRA = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(1,))
    DDEC = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(2,))
    
    coords = zip(DRA, DDEC)
    
    ra = []
    dec = []
    
    pix = wcs.wcs_world2pix(coords, 1)
    
    for i in range(0,len(pix)):
        ra.append(pix[i][0])
        dec.append(pix[i][1])
        

    fits_file = path+'fits-null.fits'
    hdu = fits.open(fits_file)[0]
    
    vcmap = vn-np.mean(vn)

    Vmin = min(vcmap)
    Vmax = max(vcmap)
    
    fig1 = plt.subplot(1, 1, 1, projection=wcs)
    plt.imshow(hdu.data, origin='lower',cmap='gray_r')
    p1=plt.scatter(ra, dec, s=30, c=vcmap, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)      
    fig1.set_ylabel('DEC')
    fig1.set_xlabel('RA')
    plt.xlim(min(ra)-0.1*max(ra), max(ra)+0.1*max(ra))
    plt.ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec))  
    cax,kw = mpl.colorbar.make_axes(fig1.axes, fraction=0.055, pad=0.07)
    cb = plt.colorbar(p1, cax=cax, **kw)
    cb.set_label('Radial Velocity (km/s)')
    plt.show() 
    
def plot_vmap_tripleplot(path, vn, cs, gals):
    RAblue = []
    DECblue= []
    RAred = []
    DECred = []
    vnb = []
    vnr = []
    
    RAall, DECall,RAb, DECb, RAr, DECr, vnbs, vnrs = [],[],[],[],[],[],[],[]
    
    for j in range(0, len(gals)):  
        
        gal = gals[j]
    
        DRA = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(1,))
        DDEC = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(2,))    
        
        col1 = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(4,))
        col2 = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(5,))
        col = col1 - col2
        
        n = len(vn[j])
        
        for i in range(0,n):
            if col[i] < cs[j]:
                RAblue.append(DRA[i])
                DECblue.append(DDEC[i])
                vnb.append(vn[j][i])
            else:
                RAred.append(DRA[i])
                DECred.append(DDEC[i])
                vnr.append(vn[j][i])
        
        RAall.append(DRA)
        DECall.append(DDEC)
        RAb.append(RAblue)
        DECb.append(DECblue)
        RAr.append(RAred)
        DECr.append(DECred)
        vnbs.append(vnb)
        vnrs.append(vnr)
        
          
    vcmaps, vrcmaps, vbcmaps, Vmins, Vmaxs, ras, decs, rars, decrs, rabs, decbs, hdus, coords_s = [],[],[],[],[],[],[],[],[],[],[],[],[]
    
    for k in range(0, len(gals)):
        
        coords = zip(RAall[k], DECall[k])
        coords_blue = zip(RAb[k], DECb[k])
        coords_red = zip(RAr[k], DECr[k])
        
        fits_file = '../Galaxies/'+gals[j]+'/'+'fits-null.fits'
        hdu = fits.open(fits_file)[0]
        wcs = WCS(hdu.header)
        
        ra = []
        dec = []
        rab = []
        decb = []
        rar = []
        decr = []
        
        pix = wcs.wcs_world2pix(coords, 1)
        pix_b = wcs.wcs_world2pix(coords_blue, 1)
        pix_r = wcs.wcs_world2pix(coords_red, 1)
        
        for i in range(0,len(pix)):
            ra.append(pix[i][0])
            dec.append(pix[i][1])
        for i in range(0,len(pix_b)):
            rab.append(pix_b[i][0])
            decb.append(pix_b[i][1])
        for i in range(0, len(pix_r)):
            rar.append(pix_r[i][0])
            decr.append(pix_r[i][1])
        
        print 'N'+gals[k]
        print 'number of blues:',len(rab)
        print 'number of reds:', len(rar)
        
        vcmap = vn[k] - np.mean(vn[k])
        vrcmap = vnrs[k] - np.mean(vnrs[k])
        vbcmap = vnbs[k] - np.mean(vnbs[k])
    
        Vmin = max(vcmap)
        Vmax = min(vcmap)
        
        vcmaps.append(vcmap)
        vrcmaps.append(vrcmap)
        vbcmaps.append(vbcmap)
        Vmins.append(Vmin)
        Vmaxs.append(Vmax)
        ras.append(ra)
        decs.append(decs)
        rabs.append(rab)
        decbs.append(decb)        
        rars.append(rars)
        decrs.append(decr)
        hdus.append(hdu.data)
        coords_s.append(coords)
    
    f, ax = plt.subplots(3, 3, projection=wcs) 
    
    for l in range(0, len(gals)):
        
        #fig1 = plt.subplot(1, 3, 1, projection=wcs)
        ax[l][1].imshow(hdu.data, origin='lower',cmap='gray_r')
        ax[l][1].scatter(ra, dec, s=30, c=vcmap, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)  
        decx = f.coords_s[l][1]
        decx.set_axislabel('DEC')
        ax[l][1].set_xlim(min(ra)-0.1*max(ra), max(ra)+0.05*max(ra))
        ax[l][1].set_ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec))
        
        #fig2 = plt.subplot(1, 3, 2, projection=wcs)
        ax[l][2].imshow(hdu.data, origin='lower',cmap='gray_r')
        ax[l][2].scatter(rar, decr, s=30, c=vrcmap, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)  
        rax = f.coords[l][0]
        rax.set_axislabel('RA')
        decx = f.coords[l][1]
        decx.set_ticklabel_visible(False)
        ax[l][2].set_xlim(min(ra)-0.1*max(ra), max(ra)+0.05*max(ra))
        ax[l][2].set_ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec))
        
        #fig3 = plt.subplot(1, 3, 3, projection=wcs)
        decx = f.coords[l][1]
        decx.set_ticklabel_visible(False)
        ax[l][3].imshow(hdu.data, origin='lower',cmap='gray_r')
        p3 = ax[l][3].scatter(rab, decb, s=30, c=vbcmap, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)   
        ax[l][3].set_xlim(min(ra)-0.1*max(ra), max(ra)+0.05*max(ra))
        ax[l][3].set_ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec)) 
        cax,kw = mpl.colorbar.make_axes(f.axes, fraction=0.047, shrink=1.0, pad=0.07)
        cb = ax[l][3].colorbar(p3, cax=cax, **kw)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        cb.set_label('Radial Velocity (km/s)')
    
        plt.subplots_adjust(wspace=.0010) 
        
    plt.savefig('vmaps_allgals.png', dpi=300, format='png')    
    plt.show() 

def find_mgc(RA, DEC, x, y):     #routine to find the Mth closest GC
    d = np.zeros(len(RA))         #array with the distance to each GC  
    
    for i in range(0, len(RA)):
        d[i] = np.sqrt((x-RA[i])**2+(y-DEC[i])**2)  # calculating distances
        
    ds = np.sort(d)     #putting distances in order to get the Mth closest
    ind = np.where(d==ds[19])  #finding the index of the Mth closest within original distance vector "d"
    
    #in this case we have M=20 (i.e. array index = 19)
    
    return RA[ind], DEC[ind]  #return RA and DEC of Mth closest GC

def aks(V, RA, DEC, x, y, A, B, dV):          #adaptive kernel smoothing (Coccatto et al, 2009)
    w = np.zeros(len(V)) #weight vector (here lies all the method!)
    
    xm, ym = find_mgc(RA, DEC, x, y)  #RA and DEC of Mth closest GC
    k = A*np.sqrt((x-xm)**2+(y-ym)**2)+B     #kernel bandwidth, based on distance to Mth GC and constants A and B  
    for i in range(0, len(V)):         
        D = np.sqrt((RA[i]-x)**2+(DEC[i]-y)**2) #given 1 GC, find the distance to each other one              
        w[i] = np.exp((-D**2)/2*k**2)            #add the ith GC weight to the weight array (so we have N weights)
    
    ss=np.zeros(len(V))      #dispersions (smoothed)
    vs=np.zeros(len(V))      #velocities (smoothed)
    
#    for i in range(0, len(V)):
#        vs[i] = ((V[i]*w[i])/sum(w))    #for each ith GC, calculate the relative veloticy
#    vs = sum(vs)
    
    vs = (sum(V*w)/sum(w))
    
    for i in range(0, len(V)):
        ss[i] = (((V[i]**2)*w[i])/(sum(w))) #calculate for ith GC relative dispersion
    ss = np.sqrt(sum(ss)-vs**2-dV**2)
      
    return vs, ss

def aks_err(V, RA, DEC, A, B, dV, vs, ss):
    n = len(RA)
    vsim = np.zeros(n)
    vsim_vs = np.array([np.zeros(n) for x in range(0, n)])
    vsim_ss = np.array([np.zeros(n) for x in range(0, n)])
    ssim = np.zeros(n)
    
    print 'calculating...'
    
    for i in range(0, n):
        rnd_ss = np.sqrt(ss[i]**2+dV**2)
        vsim[i] = vs[i]+np.random.normal(loc=0.0, scale=rnd_ss)
        
        for j in range(0, n):
            vsim_vs[i,j], vsim_ss[i,j] = aks(vsim, RA, DEC, RA[j], DEC[j], A, B, dV)
        ssim[i] = np.std(vsim_vs[i])
        
    return ssim
                 

def read_catalog_gc(galcat):
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))
    
    V = np.loadtxt(galcat, usecols=(3,))        
    
    return RA, DEC, V 
    
############    
############################################### main code 
############
    
#gal = raw_input('Enter the galaxy number: ')

galaxies = ['2768', '3115', '7457']
c_seps = []
galcats = []
paths = []
RAs, DECs, Vs = [], [], []
RAps, DECps, Vps = [], [], []
vns = []

for gal in galaxies:
    
    galinput = galaxy_inputs(gal)
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]    
    
    c_seps.append(float(inp[4]))   #colour separation
    
    galcat = '../Galaxies/'+gal+'/'+'N'+gal+'GC.dat'
    galcats.append(galcat)
    
    paths.append('../Galaxies/'+gal+'/') 

    
    RAg, DECg, V = read_catalog_gc(galcat)
    RAs.append(RAg)
    DECs.append(DECg)
    Vs.append(V)
    
    #A = 0.48
    #B = 1.53
    dV = 30.0
    vn = np.zeros(len(V))   # array with smoothed velocities
    sn = np.zeros(len(V))   # array with smoothed velocity dispersions
    pv = 0.0
    
    upt = raw_input('Estimate A and B for %s? y/n: '%gal)
    if upt == 'y':
    
        print('calculating best A and B...')
        
        it = 1
        pva = ()
        As = ()
        Bs = ()
        
        for it in range(1, 10000):
                A = np.random.uniform(0.0, 1.0)
                B = np.random.uniform(0.0, 100.0)
                for i in range(0, len(V)):        # run aks for each GC
                    vn[i], sn[i] = aks(V, RAg, DECg, RAg[i], DECg[i], A, B, dV)
                #estimating A and B:     
                KD, pv = stats.ks_2samp(V, vn)
                print it, A, B, KD, pv
                
                As = np.append(As, A)
                Bs = np.append(Bs, B)
                pva = np.append(pva,pv)
                if pv >= 0.5:
                    break
                
        print(max(pva)) 
        A = As[int(np.where(pva.max())[0])]
        B = Bs[int(np.where(pva.max())[0])]
        f = open('A_B_params'+gal, 'w')
        print >>f, A, B, max(pva)
        f.close()
    else:
        A = np.loadtxt('A_B_params'+gal, usecols=(0,))
        B = np.loadtxt('A_B_params'+gal, usecols=(1,))   
    
    
    for i in range(0, len(V)):        # run aks for each GC
        vn[i], sn[i] = aks(V, RAg, DECg, RAg[i], DECg[i], A, B, dV)    
        
    vns.append(vn)    
    
    V_c = V-np.mean(V)
    vn_c = vn-np.mean(vn)
    
    fig = plt.subplot(1, 2, 1)
    plt.text(min(RAg)+0.05,max(DECg)-0.01, 'NGC'+gal+'\n$\overline{V}_{obs}$='+str(int(np.mean(V))))
    plt.scatter(RAg, DECg, s=30, c=V_c, alpha=0.7)
    cb1 = plt.colorbar(fraction=0.04, orientation='horizontal')
    cb1.set_label('$Velocity$ $(km/s)$')
    plt.title('$V_{observed}$')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.gca().invert_xaxis()
    fig2 = plt.subplot(1, 2, 2)
    plt.scatter(RAg, DECg, s=30, c=vn_c, vmin=np.min(V_c), vmax=np.max(V_c), alpha=0.7)
    plt.text(min(RAg)+0.06,max(DECg)-0.01, '$\overline{V}_{smoothed}$='+str(int(np.mean(vn))))
    plt.title('$V_{smoothed}$')
    plt.xlabel('RA')
    cb = plt.colorbar(fraction=0.04, orientation='horizontal')
    cb.set_label('$Velocity_{smoothed}$ $(km/s)$')
    plt.gca().invert_xaxis()
    #fig2 = plt.subplot(1, 3, 3)
    #plt.scatter(RAg, DECg, s=30, c=sn, alpha=0.7)
    #plt.title('$\sigma_{smoothed}$')
    #plt.xlabel('RA')
    #cb = plt.colorbar(fraction=0.04, orientation='horizontal')
    #cb.set_label('$\sigma_{smoothed}$ $(km/s)$')
    plt.show()
    
    #if c_sep==-1000:
    #    plot_vmap_singleplot(path, vn)
    #else:
    #    plot_vmap_tripleplot(path, vn, c_sep)
    
plot_vmap_tripleplot(paths, vns, c_seps, galaxies)

ch = raw_input('run Monte Carlo sim to get errors? y/n: ')

if ch=='y':
    err = np.zeros(100)
    for i in range(0, 100):    
        std = aks_err(V, RAg, DECg, A, B, dV, vn, sn)
        err[i] = np.mean(std)
        print std        
        
    print 'standard deviation after 100 simulations: ', np.mean(err)
