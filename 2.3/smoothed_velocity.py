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
from astropy import units as u
label_size = 9
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

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
    
def plot_vmap_tripleplot(path, vn, vnp, cs, gals):
    RAblue = []
    DECblue= []
    RAred = []
    DECred = []
    vnb = []
    vnr = []

    RAall, DECall, RAb, DECb, RAr, DECr, vnbs, vnrs = [],[],[],[],[],[],[],[]
    RAp, DECp = [],[]

    for j in range(0, len(gals)):
        
        if j <= 1:
            RAblue = []
            DECblue= []
            RAred = []
            DECred = []
            vnb = []
            vnr = []
            
            gal = gals[j]
        
            DRA = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(1,))
            DDEC = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(2,))
            DRAp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(1,))
            DDECp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(2,)) 
            
            col1 = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(4,))
            col2 = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(5,))
            col = col1 - col2
            
            for i in range(0,len(DRA)):
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
            RAb.append(np.asarray(RAblue))
            DECb.append(np.asarray(DECblue))
            RAr.append(np.asarray(RAred))
            DECr.append(np.asarray(DECred))
            vnbs.append(np.asarray(vnb))
            vnrs.append(np.asarray(vnr))
            
            RAp.append(DRAp)
            DECp.append(DDECp)
            
        else:
            gal = gals[j]
        
            DRA = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(1,))
            DDEC = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'GC.dat', usecols=(2,))   
            DRAp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(1,))
            DDECp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(2,))
            
            RAall.append(DRA)
            DECall.append(DDEC)           
            RAp.append(DRAp)
            DECp.append(DDECp)
            
    RAall = np.asanyarray(RAall)
    DECall = np.asanyarray(DECall)
    RAb = np.asanyarray(RAb)
    DECb = np.asanyarray(DECb)
    RAr = np.asanyarray(RAr)
    DECr = np.asanyarray(DECr)
    vnbs = np.asanyarray(vnbs)
    vnrs = np.asanyarray(vnrs)
    RAp = np.asanyarray(RAp)
    DECp = np.asanyarray(DECp)
        
    vcmaps, vrcmaps, vbcmaps, Vmins, Vmaxs, ras, decs, rars, decrs, rabs, decbs, hdus, coords_s, wcss = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    vcmapps, raps, decps = [],[],[]   
    
    f = plt.figure()
    adj_max = [150,100]
    adj_min = [20,20]
    for k in range(0, len(gals)):
        if k <= 1:
        
            coords = zip(RAall[k], DECall[k])
            coords_blue = zip(RAb[k], DECb[k])
            coords_red = zip(RAr[k], DECr[k])
            coords_p = zip(RAp[k], DECp[k])
            
            fits_file = '../Galaxies/'+gals[k]+'/'+'fits-null.fits'
            hdu = fits.open(fits_file)[0]
            wcs = WCS(hdu.header)
            
            ra = []
            dec = []
            rab = []
            decb = []
            rar = []
            decr = []
            rap = []
            decp = []
            
            pix = wcs.wcs_world2pix(coords, 1)
            pix_p = wcs.wcs_world2pix(coords_p, 1)
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
            for i in range(0, len(pix_p)):
                rap.append(pix_p[i][0])
                decp.append(pix_p[i][1])
            
            vcmap = vn[k] - np.mean(vn[k])
            vrcmap = vnrs[k] - np.mean(vnrs[k])
            vbcmap = vnbs[k] - np.mean(vnbs[k])
            vcmapp = vnp[k] - np.mean(vnp[k])
        
            Vmin = min(vcmap)+adj_min[k]
            Vmax = max(vcmap)-adj_max[k]
            
            print 'N'+gals[k], max(vcmap)-np.median(vcmap), np.median(vcmap)-min(vcmap)
            
            vcmapps.append(vcmapp)
            vcmaps.append(vcmap)
            vrcmaps.append(vrcmap)
            vbcmaps.append(vbcmap)
            Vmins.append(Vmin)
            Vmaxs.append(Vmax)
            raps.append(rap)
            decps.append(decp)
            ras.append(ra)
            decs.append(dec)
            rabs.append(rab)
            decbs.append(decb)        
            rars.append(rar)
            decrs.append(decr)
            hdus.append(hdu.data)
            coords_s.append(coords)
            wcss.append(wcs)   
    
            fig1 = f.add_subplot(3, 4, 4*k+1, projection=wcss[k])
            fig1.imshow(hdus[k], origin='lower',cmap='gray_r')
            fig1.text(min(ras[k])-0.08*max(ras[k]),max(decs[k]), 'All GCs', fontsize=10)
            fig1.scatter(ras[k], decs[k], s=15, c=vcmaps[k], vmin=Vmins[k], vmax=Vmaxs[k], edgecolors='None', alpha=1.0)  
            #fig1.scatter(ras[l], decs[l], s=30, c=vcmaps[l], edgecolors='None', alpha=1.0)  
            fig1.set_xlim(min(ras[k])-0.1*max(ras[k]), max(ras[k])+0.05*max(ras[k]))
            fig1.set_ylim(min(decs[k])-0.1*max(decs[k]), max(decs[k])+0.1*max(decs[k]))
            rax1 = fig1.coords[0]  
            rax1.set_axislabel(' ', minpad=-10.0)
            rax1.set_major_formatter('dd:mm')
            rax1.set_separator(('d', "'"))
            #rax1.set_ticks(spacing=6.*u.arcmin)
            rax1.set_ticks(number=2)
            decx1 = fig1.coords[1]
            decx1.set_axislabel('DEC', fontsize=15, minpad=0.2)
            decx1.set_major_formatter('dd:mm')
            decx1.set_separator(('d', "'"))
            decx1.set_ticks(number=3)
            rax1.display_minor_ticks(True)
            decx1.display_minor_ticks(True)
            rax1.set_ticklabel(size=9)
            decx1.set_ticklabel(size=9)
            
            fig2 = f.add_subplot(3, 4, 4*k+2, projection=wcss[k])
            fig2.imshow(hdus[k], origin='lower',cmap='gray_r')
            fig2.text(min(ras[k])-0.08*max(ras[k]),max(decs[k]), 'Red GCs', fontsize=10)
            fig2.scatter(rars[k], decrs[k], s=15, c=vrcmaps[k], vmin=Vmins[k], vmax=Vmaxs[k], edgecolors='None', alpha=1.0)  
            #fig2.scatter(rars[l], decrs[l], s=30, c=vrcmaps[l], edgecolors='None', alpha=1.0)  
            fig2.set_xlim(min(ras[k])-0.1*max(ras[k]), max(ras[k])+0.05*max(ras[k]))
            fig2.set_ylim(min(decs[k])-0.1*max(decs[k]), max(decs[k])+0.1*max(decs[k]))
            rax2 = fig2.coords[0]
            rax2.set_major_formatter('dd:mm')
            rax2.set_separator(('d', "'"))
            rax2.set_ticks(number=2)
            #rax2.set_ticks(spacing=6.*u.arcmin)
            decx2 = fig2.coords[1]
            decx2.set_ticklabel_visible(False)
            rax2.display_minor_ticks(True)
            decx2.display_minor_ticks(True)
            rax2.set_ticklabel(size=9)
            decx2.set_ticklabel(size=9)

            
            fig3 = f.add_subplot(3, 4, 4*k+3, projection=wcss[k])
            fig3.imshow(hdus[k], origin='lower',cmap='gray_r')
            fig3.text(min(ras[k])-0.08*max(ras[k]),max(decs[k]), 'Blue GCs', fontsize=10)            
            fig3.scatter(rabs[k], decbs[k], s=15, c=vbcmaps[k], vmin=Vmins[k], vmax=Vmaxs[k], edgecolors='None', alpha=1.0)
            #p3 = fig3.scatter(rabs[l], decbs[l], s=30, c=vbcmaps[l], edgecolors='None', alpha=1.0)    
            #cax,kw = mpl.colorbar.make_axes(fig3.axes, fraction=0.047, shrink=1.0, pad=0.07)
            fig3.set_xlim(min(ras[k])-0.1*max(ras[k]), max(ras[k])+0.05*max(ras[k]))
            fig3.set_ylim(min(decs[k])-0.1*max(decs[k]), max(decs[k])+0.1*max(decs[k]))
            #fig3.text(max(ras[k])-0.2*max(ras[k]),max(decs[k]), 'N'+gals[k], fontsize=12)
            rax3 = fig3.coords[0]      
            rax3.set_axislabel(' ', minpad=-30.0)
            rax3.set_major_formatter('dd:mm')
            rax3.set_separator(('d', "'"))
            #rax3.set_ticks(spacing=6.*u.arcmin)
            rax3.set_ticks(number=2)
            decx3 = fig3.coords[1]
            decx3.set_ticklabel_visible(False)
            rax3.display_minor_ticks(True)
            decx3.display_minor_ticks(True)
            rax3.set_ticklabel(size=9)
            decx3.set_ticklabel(size=9)
              
            fig4 = plt.subplot(3, 4, 4*k+4, projection=wcs)
            fig4.imshow(hdu.data, origin='lower',cmap='gray_r')
            fig4.text(min(ras[k])-0.08*max(ras[k]),max(decs[k]), 'PNe', fontsize=10)
            fig4.text(max(ras[k])-0.30*max(ras[k]),max(decs[k]), 'N'+gals[k], fontsize=10)            
            p4 = fig4.scatter(raps[k], decps[k], s=15, c=vcmapps[k], vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)      
            fig4.set_xlim(min(ras[k])-0.1*max(ras[k]), max(ras[k])+0.05*max(ras[k]))
            fig4.set_ylim(min(decs[k])-0.1*max(decs[k]), max(decs[k])+0.1*max(decs[k]))
            rax4 = fig4.coords[0]  
            rax4.set_axislabel(' ', minpad=-10.0)
            rax4.set_major_formatter('dd:mm')
            rax4.set_separator(('d', "'"))
            #rax1.set_ticks(spacing=6.*u.arcmin)
            rax4.set_ticks(number=2)
            decx4 = fig4.coords[1]
            decx4.set_ticklabel_visible(False)
            rax4.display_minor_ticks(True)
            decx4.display_minor_ticks(True) 
            #rax4.set_axislabel('RA', fontsize=15, minpad=0.4)    
            rax4.set_ticklabel(size=9)
            decx4.set_ticklabel(size=9)
            if k==0:
                cbaxes = f.add_axes([0.9, 0.645, 0.02, 0.254])
                cb = plt.colorbar(p4, cax = cbaxes, orientation='vertical')
                #cb = plt.colorbar(p4, orientation='vertical', fraction=0.07)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cb.locator = tick_locator
                cb.update_ticks()
                cb.set_label('Velocity (km/s)', fontsize=10)
            if k==1:
                cbaxes = f.add_axes([0.9, 0.391, 0.02, 0.218])
                cb = plt.colorbar(p4, cax = cbaxes, orientation='vertical')
                #cb = plt.colorbar(p4, orientation='vertical', fraction=0.07)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cb.locator = tick_locator
                cb.update_ticks()
                cb.set_label('Velocity (km/s)', fontsize=10)
            #cbaxes.xaxis.set_ticks_position('top')
            #cbaxes.xaxis.set_label_position('top')
            
            if k==1:
                #rax2.set_axislabel('RA', fontsize=15, minpad=0.4)
                rax3.set_axislabel('RA', fontsize=15, minpad=0.4)
                rax4.set_axislabel('RA', fontsize=15, minpad=0.4)
        else:
            hdu_list = fits.open(path[k]+'fits-null.fits')
            hdu = hdu_list[0]
            wcs = WCS(hdu.header)
            
            DRA = np.loadtxt('../Galaxies/'+gals[k]+'/'+'N'+gals[k]+'GC.dat', usecols=(1,))
            DDEC = np.loadtxt('../Galaxies/'+gals[k]+'/'+'N'+gals[k]+'GC.dat', usecols=(2,))
            DRAp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(1,))
            DDECp = np.loadtxt('../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat', usecols=(2,))
            
            coords = zip(DRA, DDEC)
            coords_p = zip(DRAp, DDECp)
            
            ra = []
            dec = []
            
            pix = wcs.wcs_world2pix(coords, 1)
            
            rap = []
            decp = []
            
            pix_p = wcs.wcs_world2pix(coords_p, 1)
            
            for i in range(0,len(pix)):
                ra.append(pix[i][0])
                dec.append(pix[i][1])
            for i in range(0, len(pix_p)):
                rap.append(pix_p[i][0])
                decp.append(pix_p[i][1])                
        
            fits_file = path[k]+'fits-null.fits'
            hdu = fits.open(fits_file)[0]
            
            vcmap = vn[k]-np.mean(vn[k])
            vcmapp = vnp[k]-np.mean(vnp[k])
        
            Vmin = min(vcmap)
            Vmax = max(vcmap)
            
            print 'N'+gals[k], max(vcmap)-np.median(vcmap), np.median(vcmap)-min(vcmap)
            
            fig1 = plt.subplot(3, 4, 9, projection=wcs)
            fig1.imshow(hdu.data, origin='lower',cmap='gray_r')
            fig1.text(min(ra)-0.08*max(ra),max(dec)+0.05*max(dec), 'All GCs', fontsize=10)
            fig1.scatter(ra, dec, s=15, c=vcmap, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)      
            fig1.set_xlim(min(ra)-0.1*max(ra), max(ra)+0.1*max(ra))
            fig1.set_ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec))
            rax = fig1.coords[0]  
            rax.set_axislabel(' ', minpad=-10.0)
            rax.set_major_formatter('dd:mm')
            rax.set_separator(('d', "'"))
            #rax1.set_ticks(spacing=6.*u.arcmin)
            rax.set_ticks(number=3)
            decx = fig1.coords[1]
            decx.set_axislabel('DEC', fontsize=15, minpad=0.2)
            decx.set_major_formatter('dd:mm')
            decx.set_separator(('d', "'"))
            decx.set_ticks(number=3)
            rax.display_minor_ticks(True)
            decx.display_minor_ticks(True) 
            rax.set_axislabel('RA', fontsize=15, minpad=0.4)
            rax.set_ticklabel(size=9)
            decx.set_ticklabel(size=9)

            fig4 = plt.subplot(3, 4, 10, projection=wcs)
            fig4.imshow(hdu.data, origin='lower',cmap='gray_r')
            fig4.text(min(ra)-0.08*max(ra),max(dec)+0.05*max(dec), 'PNe', fontsize=10)
            fig4.text(max(ra)-0.1*max(ra), max(dec)+0.05*max(dec), 'N'+gals[k], fontsize=10)
            p4 = fig4.scatter(rap, decp, s=15, c=vcmapp, vmin=Vmin, vmax=Vmax, edgecolors='None', alpha=1.0)   
            fig4.set_xlim(min(ra)-0.1*max(ra), max(ra)+0.1*max(ra))
            fig4.set_ylim(min(dec)-0.1*max(dec), max(dec)+0.1*max(dec))
            rax4 = fig4.coords[0]  
            rax4.set_axislabel(' ', minpad=-10.0)
            rax4.set_major_formatter('dd:mm')
            rax4.set_separator(('d', "'"))
            #rax1.set_ticks(spacing=6.*u.arcmin)
            rax4.set_ticks(number=3)
            decx4 = fig4.coords[1]
            decx4.set_ticklabel_visible(False)
            rax4.display_minor_ticks(True)
            decx4.display_minor_ticks(True) 
            rax4.set_axislabel('RA', fontsize=15, minpad=0.4)     
            rax4.set_ticklabel(size=9)
            decx4.set_ticklabel(size=9)
            cbaxes = f.add_axes([0.511, 0.12, 0.02, 0.214])
            cb = plt.colorbar(p4, cax = cbaxes, orientation='vertical')
            #cb = plt.colorbar(p4, orientation='vertical', fraction=0.07)
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            cb.update_ticks()
            cb.set_label('Velocity (km/s)', fontsize=10)
            
    f.subplots_adjust(wspace=.0)
    f.subplots_adjust(hspace=.07)        
        
    plt.savefig('vmaps_allgals.png', dpi=300, format='png')    
    plt.show() 

def find_mgc(RA, DEC, x, y):     #routine to find the Mth closest GC
    d = np.zeros(len(RA))         #array with the distance to each GC  
    
    for i in range(0, len(RA)):
        d[i] = np.sqrt((x-RA[i])**2+(y-DEC[i])**2)  # calculating distances
        
    ds = np.sort(d)     #putting distances in order to get the Mth closest
    ind = np.where(d==ds[19])  #finding the index of the Mth closest within original distance vector "d"
    if len(ind[0])>1:
        ind=ind[0][0]
    
    #in this case we have M=20 (i.e. array index = 19)
    
    return RA[ind], DEC[ind]  #return RA and DEC of Mth closest GC

def aks(Vas, RA, DEC, x, y, A, B, dV):          #adaptive kernel smoothing (Coccatto et al, 2009)
    w = np.zeros(len(Vas)) #weight vector (here lies all the method!)
    
    xm, ym = find_mgc(RA, DEC, x, y)  #RA and DEC of Mth closest GC
    k = A*np.sqrt((x-xm)**2+(y-ym)**2)+B     #kernel bandwidth, based on distance to Mth GC and constants A and B  
    for i in range(0, len(Vas)):         
        D = np.sqrt((RA[i]-x)**2+(DEC[i]-y)**2) #given 1 GC, find the distance to each other one             
        w[i] = np.exp((-D**2)/2*k**2)            #add the ith GC weight to the weight array (so we have N weights)
    
    ss=np.zeros(len(Vas))      #dispersions (smoothed)
    vs=np.zeros(len(Vas))      #velocities (smoothed)
    
#    for i in range(0, len(V)):
#        vs[i] = ((V[i]*w[i])/sum(w))    #for each ith GC, calculate the relative veloticy
#    vs = sum(vs)
    
    vs = (sum(Vas*w)/sum(w))
    
    for i in range(0, len(Vas)):
        ss[i] = (((Vas[i]**2)*w[i])/(sum(w))) #calculate for ith GC relative dispersion
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
galcatsp = []
paths = []
RAs, DECs, Vs = [], [], []
RAps, DECps, Vps = [], [], []
vns = []
vnps = []

for gal in galaxies:
    
    galinput = galaxy_inputs(gal)
    
    with open(galinput, 'r') as f:
        inp = [x.split(' ')[0] for x in open(galinput).readlines()]    
    
    c_seps.append(float(inp[4]))   #colour separation
    
    galcat = '../Galaxies/'+gal+'/'+'N'+gal+'GC.dat'
    galcatp = '../Galaxies/'+gal+'/'+'N'+gal+'PNE.dat'
    galcats.append(galcat)
    galcatsp.append(galcatp)
    
    paths.append('../Galaxies/'+gal+'/') 

    
    RAg, DECg, V = read_catalog_gc(galcat)
    RAp, DECp, Vp = read_catalog_gc(galcatp)
    
    RAs.append(RAg)
    DECs.append(DECg)
    Vs.append(V)
    RAps.append(RAp)
    DECps.append(DECp)
    Vps.append(Vp)
    
    #A = 0.48
    #B = 1.53
    dV = 30.0
    vn = np.zeros(len(V))   # array with smoothed velocities
    vnp = np.zeros(len(Vp))
    sn = np.zeros(len(V))   # array with smoothed velocity dispersions
    snp = np.zeros(len(Vp))
    pv = 0.0
    
    #upt = raw_input('Estimate A and B for %s? y/n: '%gal)
    upt='n'
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
    for i in range(0, len(Vp)):        # run aks for each GC
        vnp[i], snp[i] = aks(Vp, RAp, DECp, RAp[i], DECp[i], A, B, dV)    

#    vnp = Vp
        
    vns.append(vn)   
    vnps.append(vnp) 
    
#    V_c = V-np.mean(V)
#    vn_c = vn-np.mean(vn)
    
#    fig = plt.subplot(1, 2, 1)
#    plt.text(min(RAg)+0.05,max(DECg)-0.01, 'NGC'+gal+'\n$\overline{V}_{obs}$='+str(int(np.mean(V))))
#    plt.scatter(RAg, DECg, s=30, c=V_c, alpha=0.7)
#    cb1 = plt.colorbar(fraction=0.04, orientation='horizontal')
#    cb1.set_label('$Velocity$ $(km/s)$')
#    plt.title('$V_{observed}$')
#    plt.xlabel('RA')
#    plt.ylabel('DEC')
#    plt.gca().invert_xaxis()
#    fig2 = plt.subplot(1, 2, 2)
#    plt.scatter(RAg, DECg, s=30, c=vn_c, vmin=np.min(V_c), vmax=np.max(V_c), alpha=0.7)
#    plt.text(min(RAg)+0.06,max(DECg)-0.01, '$\overline{V}_{smoothed}$='+str(int(np.mean(vn))))
#    plt.title('$V_{smoothed}$')
#    plt.xlabel('RA')
#    cb = plt.colorbar(fraction=0.04, orientation='horizontal')
#    cb.set_label('$Velocity_{smoothed}$ $(km/s)$')
#    plt.gca().invert_xaxis()
#    #fig2 = plt.subplot(1, 3, 3)
#    #plt.scatter(RAg, DECg, s=30, c=sn, alpha=0.7)
#    #plt.title('$\sigma_{smoothed}$')
#    #plt.xlabel('RA')
#    #cb = plt.colorbar(fraction=0.04, orientation='horizontal')
#    #cb.set_label('$\sigma_{smoothed}$ $(km/s)$')
#    plt.show()
    
    #if c_sep==-1000:
    #    plot_vmap_singleplot(path, vn)
    #else:
    #    plot_vmap_tripleplot(path, vn, c_sep)
    
plot_vmap_tripleplot(paths, vns, vnps, c_seps, galaxies)

#ch = raw_input('run Monte Carlo sim to get errors? y/n: ')
ch='n'
if ch=='y':
    err = np.zeros(100)
    for i in range(0, 100):    
        std = aks_err(V, RAg, DECg, A, B, dV, vn, sn)
        err[i] = np.mean(std)
        print std        
        
    print 'standard deviation after 100 simulations: ', np.mean(err)
