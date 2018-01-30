import numpy as np
import pylab
import uncertainties.unumpy as unp
import matplotlib.ticker as plticker
from bin_functions import galaxy_inputs

#TO DO: This code is too big and repetitive. We should create another function code to import, this time for plotting things.

def read_catalog(cat_name):
    GCdens = np.loadtxt(cat_name, usecols=(2,))
    GCdenserr = np.loadtxt(cat_name, usecols=(3,))
    rgc = np.loadtxt(cat_name, usecols=(6,))
    binsize_gc = np.loadtxt(cat_name, usecols=(8,))
    
    return GCdens, GCdenserr, rgc, binsize_gc
    
def read_catalog_log(cat_name):
    N = np.loadtxt(cat_name, usecols=(5,))
    area = np.loadtxt(cat_name, usecols=(4,))
    r = np.loadtxt(cat_name, usecols=(6,))
    
    return N, area, r
    
def log_galfit(x, mue, mud, bn, re, n, h):
    re_min = re/60
    hmin = h/60
    xmin = x/60
    srsc2 = mue*np.exp(-bn*((xmin/re_min)**(1/n)-1))
    exp12 = mud*np.exp(-xmin/hmin)
    total2 = srsc2+exp12   
 
    xlog = np.log10(xmin)
    srsclog = np.log10(srsc2)
    exp1log = np.log10(exp12)
    totallog = np.log10(total2)

    return xlog, srsclog, exp1log, totallog    
    
def galfit(x, mue, mud, bn, re, n, h):    

    srsc = mue*np.exp(-bn*((x/re)**(1/n)-1))
    exp1 = mud*np.exp(-x/h)
    total = srsc+exp1
        
    srsc = np.log10(srsc)
    exp1 = np.log10(exp1)
    total = np.log10(total)
    
    return srsc, exp1, total
    
    
def config_logplots(N, area, r):
    area = area/3600
    
    density = N/area
    denslog = np.log10(density)
    
    err = np.sqrt(N)/area
    errlog = err/(2.3*(density))
    
    rmin = r/60
        
    rlog = np.log10(rmin)
    #TO DO: add binsizes

    return denslog, errlog, rlog   
    
def read_photometric_catalog(cat_name):
    phot_rad = np.loadtxt(cat_name, usecols=(0,))
    phot_gc_dens = np.loadtxt(cat_name, usecols=(1,))
    
    return phot_gc_dens, phot_rad
    
def phot_error(cat_name, phot_gc_dens, ch):
    phot_err = np.loadtxt(cat_name, usecols=(2,))
    dens_err = unp.uarray(phot_gc_dens, phot_err)
    
    if ch == 'l':    
        denslog_err = unp.log10(dens_err)
    else: 
        denslog_err = -2.5*unp.log10(dens_err)
    error = unp.std_devs(denslog_err)
    
    return error

with open('current_galaxy.dat', 'r') as f:     
    gal = [line.split()[0] for line in f]
    gal = gal[0]

galput = galaxy_inputs(gal)
inp = [x.split(' ')[0] for x in open(galput).readlines()]
    
c_sep = float(inp[4])    
    
n = float(inp[5])
bn = 2.0*n-0.327
re = float(inp[6])
h = float(inp[7])
mue = float(inp[8])
mud = float(inp[9])

#my densities of GC's

GCdensfile = 'test_GCdensity-circularbins.dat'
redGCdensfile = 'test_GCdensity-REDS.dat'
blueGCdensfile = 'test_GCdensity-BLUES.dat'

GCdens, GCdenserr, rgc, binsize_gc = read_catalog(GCdensfile)

if c_sep != -1000:
    redGCdens, redGCdenserr, rgc_reds, binsize_gc_reds = read_catalog(redGCdensfile)
    blueGCdens, blueGCdenserr, rgc_blues, binsize_gc_blues = read_catalog(blueGCdensfile)

#setting an array to plot light profiles
x = np.arange(0,1500,5.17)
nr = x.size

#plots

op = raw_input('Do we have Pota photometric GCs for 2768? y/n: ')
op2 = raw_input('Do we have PNe? y/n: ')
ch = raw_input('Type l for arcmin or r for arcsec: ')
sc = raw_input('Do you want to change the y scale of the plots? y/n: ')
if sc == 'y':    
    scale = raw_input('Enter the amount for GC: ')
    scale = float(scale)
    spne = raw_input('Enter the amount for PNe: ')
    spne = float(spne)
    
else:
    scale = 0.0
    spne = 0.0

if op == 'y':    
    photblue = gal+'/'+'SD_blue.dat'
    photred = gal+'/'+'SD_red.corr.dat'    
    
    phot_gc_dens_r, phot_rad_r = read_photometric_catalog(photred)
    phot_gc_dens_b, phot_rad_b = read_photometric_catalog(photblue)    
    
    #######log plots
    if ch == 'l':    
        N, area, rall = read_catalog_log(GCdensfile)
        N_r, area_r, rreds = read_catalog_log(redGCdensfile)
        N_b, area_b, rblues = read_catalog_log(blueGCdensfile)       
        denslog, errlog, rlog = config_logplots(N, area, rall)    
        denslog_r, errlog_r, rlog_r = config_logplots(N_r, area_r, rreds)
        denslog_b, errlog_b, rlog_b = config_logplots(N_b, area_b, rblues)
        
        phot_err_r = phot_error(photred, phot_gc_dens_r, ch)
        phot_err_b = phot_error(photblue, phot_gc_dens_b, ch)    
        phot_rad_r = np.log10(phot_rad_r)
        phot_rad_b = np.log10(phot_rad_b)        
        phot_gc_dens_r = np.log10(phot_gc_dens_r)
        phot_gc_dens_b = np.log10(phot_gc_dens_b)
        
        xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)    
       
        if op2 == 'y':
            
            
            #PNe 

            PNedensfile = 'PNEdensity-circularbins.dat'
            
            dens, denserr, rmed, binsize = read_catalog(PNedensfile)

            Npne, areapne, rpne = read_catalog_log(PNedensfile)
            denslogpne, errlogpne, rlogpne = config_logplots(Npne, areapne, rpne)    
            
            plt_all = pylab.subplot(1,1,1)
            pylab.plot(xlog, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
            pylab.plot(xlog, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
            pylab.plot(xlog, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
            pylab.errorbar(rlogpne, denslogpne+spne, yerr=errlogpne, capsize=6, capthick=1, marker = 'o',color='LightSeaGreen', markersize=8, label='PNe Density', linestyle='none')
            pylab.errorbar(rlog, denslog+scale, yerr=errlog, capsize=6, capthick=1, marker = 'o', color='magenta', markersize=10, label='All GC', linestyle='none')
            pylab.errorbar(rlog_r, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=1, marker = 'o', color='red', markersize=10, label='Red GC', linestyle='none')
            pylab.errorbar(rlog_b, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=1,marker = 'o', color='blue', markersize=10, label='Blue GC', linestyle='none')        
            pylab.legend(loc='lower left', prop={'size':10})    
            #pylab.xlim(-1,1.5)
            #pylab.ylim(-4,2.5)
            pylab.xlabel('log(R)($arcmin$)')
            pylab.ylabel('log($\Sigma$)($N/arcmin^2$)')    
            pylab.show()
        
    ################################################################################################################    
        
        plt = pylab.subplot(3,1,1)
        pylab.plot(xlog, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
        pylab.plot(xlog, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
        pylab.plot(xlog, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
        pylab.errorbar(rlog, denslog+scale, yerr=errlog, capsize=6, capthick=1, ecolor='magenta', marker = 'o', markeredgecolor='magenta', markersize=10, mfc='None', label='All GC', linestyle='none')
        pylab.legend(loc='lower left', prop={'size':10})
        #pylab.xlim(-1,1.5)
        #pylab.ylim(-4,2.5)
        
        plt2 = pylab.subplot(3,1,2)    
        pylab.plot(xlog, srsclog, color='red', linestyle='--')
        pylab.plot(xlog, exp1log, color ='blue', linestyle='--')
        pylab.plot(xlog, totallog, color='gray', linestyle='-')
        pylab.errorbar(rlog_r, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=1, ecolor='red', marker = 'o', markeredgecolor='red', markersize=10, mfc='None', label='Red GC', linestyle='none')
        pylab.errorbar(phot_rad_r, phot_gc_dens_r+scale,  yerr=phot_err_r, capsize=6, capthick=1, marker='D', color='orange', markersize=8, label='Photometric GC density (Reds)', mfc='None', linestyle= 'none')    
        pylab.legend(loc='lower left', prop={'size':10})    
        #pylab.xlim(-1,1)
        #pylab.ylim(-4,2.5)
        
        plt3 = pylab.subplot(3,1,3)
        pylab.plot(xlog, srsclog, color='red', linestyle='--')
        pylab.plot(xlog, exp1log, color ='blue', linestyle='--')
        pylab.plot(xlog, totallog, color='gray', linestyle='-')
        pylab.errorbar(rlog_b, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=1, ecolor='blue',marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='Blue GC', linestyle='none')    
        pylab.errorbar(phot_rad_b, phot_gc_dens_b+scale, yerr=phot_err_b, capsize=6, capthick=1, marker='D', color='orange', markersize=8, label='Photometric GC density (Blues)', mfc='None', linestyle= 'none')    
        pylab.legend(loc='lower left', prop={'size':10})    
        #pylab.xlim(-1,1)
        #pylab.ylim(-4,2.5)
        
        pylab.subplots_adjust(hspace=.0001)
        pylab.setp(plt.get_xticklabels(), visible=False)
        pylab.setp(plt2.get_xticklabels(), visible=False)
        plt2.set_ylabel('log($\Sigma$) ($N/arcmin^2$)', fontsize=14)
        plt3.set_xlabel('log(R) ($arcmin$)', fontsize=14)
        pylab.show()
                    
    else:
        #### regular plots
    
        #### photometric GCs ###   
        
        phot_rad_b = phot_rad_b*60
        phot_gc_dens_b = phot_gc_dens_b/3600
        
        phot_rad_r = phot_rad_r*60
        phot_gc_dens_r = phot_gc_dens_r/3600
        
        phot_err_r = phot_error(photred, phot_gc_dens_r)
        phot_err_b = phot_error(photblue, phot_gc_dens_b)    
        phot_err_r = phot_err_r/3600   
        phot_err_b = phot_err_b/3600    
        
        phot_gc_dens_r = -2.5*np.log10(phot_gc_dens_r)
        phot_gc_dens_b = -2.5*np.log10(phot_gc_dens_b)
    
        #scale factors
        
        dens = dens-8.0
        GCdens = GCdens-10.0
        redGCdens = redGCdens-10.0
        blueGCdens = blueGCdens-10.0
        phot_gc_dens_b = phot_gc_dens_b-10.0
        phot_gc_dens_r = phot_gc_dens_r-10.0
        
        srsc, exp1, total = galfit(x, mue, mud, bn, re, n, h)    
        
        #pylab.plot(rad, logrr, marker='D', linestyle='none', markersize=5, color='red', label='Ellipse Fit', mfc='None')
        plt = pylab.subplot(3,1,1)
        pylab.plot(x, srsc, color='green', label='Sersic Profile')
        pylab.plot(x, exp1, color ='c', label='Exp Disk Profile')
        pylab.plot(x, total, color='black', linestyle='--', label = 'Total Light Profile')
        pylab.errorbar(rgc, GCdens, xerr=binsize_gc, yerr=GCdenserr, marker = 'o',color='purple', markersize=10, mfc='None', label='total GC density', linestyle='none')
        pylab.legend(loc='upper right', prop={'size':10})
        #pylab.xlim(0,1100)
        #pylab.ylim(-10,10)
        
        
        plt2 = pylab.subplot(3,1,2, sharex=plt, sharey=plt)
        pylab.plot(x, srsc, color='green')
        pylab.plot(x, exp1, color ='c')
        pylab.plot(x, total, color='black', linestyle='--')
        pylab.errorbar(rgc_blues, blueGCdens, xerr=binsize_gc_blues,yerr=blueGCdenserr, marker = 'o',color='blue', markersize=10, label='blue GC Density', mfc='None', linestyle='none')
        pylab.errorbar(phot_rad_b, phot_gc_dens_b, yerr=phot_err_b, marker='D', color='orange', markersize=8, label='Photometric GC density (Blues)', mfc='None', linestyle= 'none')
        pylab.legend(loc='upper right', prop={'size':10})
        #pylab.xlim(0,1100)
        #pylab.ylim(-10,10)
        
        plt3 = pylab.subplot(3,1,3, sharex=plt, sharey=plt)
        pylab.plot(x, srsc, color='green')
        pylab.plot(x, exp1, color ='c')
        pylab.plot(x, total, color='black', linestyle='--')
        pylab.errorbar(rgc_reds, redGCdens, xerr=binsize_gc_reds, yerr=redGCdenserr,  marker = 'o',color='red', markersize=10, label='red GC Density', mfc='None', linestyle= 'none')
        pylab.errorbar(phot_rad_r, phot_gc_dens_r,  yerr=phot_err_r, marker='D', color='orange', markersize=8, label='Photometric GC density (Reds)', mfc='None', linestyle= 'none')
        pylab.legend(loc='upper right', prop={'size':10})
        #pylab.xlim(0,1100)
       # pylab.ylim(-10,10)
        pylab.plt.gca().invert_yaxis()
        
        #plt4 = pylab.subplot(4,1,4, sharex=plt, sharey=plt)
        #pylab.plot(x, srsc, color='green')
        #pylab.plot(x, exp1, color ='c')
        #pylab.plot(x, total, color='black', linestyle='--')
        #pylab.errorbar(rmed, dens, yerr=denserr, xerr=binsize, marker = 'o',color='green', markersize=8, mfc='None', label='PNe Density', linestyle='none')
        #pylab.legend(loc='upper right', prop={'size':10})
        #pylab.xlim(0,800)
        #pylab.ylim(-10,10)
        #pylab.subplots_adjust(hspace=.0001)
        #pylab.plt.gca().invert_yaxis()
        
        loc = plticker.MultipleLocator(base = 4.0)
        plt.yaxis.set_major_locator(loc)
        pylab.subplots_adjust(hspace=.0001)
        pylab.setp(plt.get_xticklabels(), visible=False)
        pylab.setp(plt2.get_xticklabels(), visible=False)
        #pylab.setp(plt3.get_xticklabels(), visible=False)
        plt2.set_ylabel('Surface Brightness ($mag/arcsec^2$)', fontsize=14)
        plt3.set_xlabel('R ($arcsec$)', fontsize=14)
        #plt.set_title('NGC2768')
        pylab.show()
else:
        #######log plots      
    
    if ch == 'l':  
        if c_sep == -1000:            #for the case where there is no bimodality
            N, area, rall = read_catalog_log(GCdensfile)     
            denslog, errlog, rlog = config_logplots(N, area, rall)    
            
            xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)    
            
            if op2 == 'y':
                            #PNe 

                PNedensfile = 'PNEdensity-circularbins.dat'
            
                dens, denserr, rmed, binsize = read_catalog(PNedensfile)
                
                Npne, areapne, rpne = read_catalog_log(PNedensfile)
                denslogpne, errlogpne, rlogpne = config_logplots(Npne, areapne, rpne) 
                
                plt = pylab.subplot(1,1,1)
                pylab.plot(x/60, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
                pylab.plot(x/60, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
                pylab.plot(x/60, totallog, color='gray', linestyle='-', label = 'Light Profile from GALFIT')
                pylab.errorbar(rpne/60, denslogpne+scale, yerr=errlogpne, capsize=6, capthick=1, ecolor='green', marker = 'o', markeredgecolor='green', markersize=10, mfc='None', label='GMOS Sample', linestyle='none')
                pylab.errorbar(rall/60, denslog+scale, yerr=errlog, capsize=6, capthick=1, ecolor='blue', marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='WFC3 Sample', linestyle='none')                
                pylab.legend(loc='lower left', prop={'size':10})
                #pylab.xlim(-1,1)
                #pylab.ylim(-4,2.5)
                plt.set_ylabel('log($\Sigma$) ($N/arcmin^2$)', fontsize=14)
                plt.set_xlabel('log(R) ($arcmin$)', fontsize=14)
                pylab.show()
            
            else:     
                plt = pylab.subplot(1,1,1)
                pylab.plot(x/60, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
                pylab.plot(x/60, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
                pylab.plot(x/60, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
                pylab.errorbar(rall/60, denslog+scale, yerr=errlog, capsize=6, capthick=1, ecolor='blue', marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='GC', linestyle='none')
                pylab.legend(loc='lower left', prop={'size':10})
                #pylab.xlim(-1,1)
                #pylab.ylim(-4,2.5)
                plt.set_ylabel('log($\Sigma$) ($N/arcmin^2$)', fontsize=14)
                plt.set_xlabel('R ($arcmin$)', fontsize=14)
                pylab.show()
            
        else:
            N, area, rall = read_catalog_log(GCdensfile)
            N_r, area_r, rreds = read_catalog_log(redGCdensfile)
            N_b, area_b, rblues = read_catalog_log(blueGCdensfile)       
            denslog, errlog, rlog = config_logplots(N, area, rall)    
            denslog_r, errlog_r, rlog_r = config_logplots(N_r, area_r, rreds)
            denslog_b, errlog_b, rlog_b = config_logplots(N_b, area_b, rblues)
            
            xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)    
            
            if op2 == 'y':
                #optional: PNE
                #PNe 

                PNedensfile = 'PNEdensity-circularbins.dat'
            
                dens, denserr, rmed, binsize = read_catalog(PNedensfile)            
                Npne, areapne, rpne = read_catalog_log(PNedensfile)
                denslogpne, errlogpne, rlogpne = config_logplots(Npne, areapne, rpne)    
                
                plt_all = pylab.subplot(1,1,1)
                pylab.plot(x/60, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
                pylab.plot(x/60, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
                pylab.plot(x/60, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
                pylab.errorbar(rpne/60, denslogpne+spne, yerr=errlogpne, capsize=6, capthick=1, marker = 'o',color='LightSeaGreen', markersize=8, label='PNe Density', linestyle='none')
                #pylab.errorbar(rlog, denslog+scale, yerr=errlog, capsize=6, capthick=1, marker = 'o', color='magenta', markersize=10, label='All GC', linestyle='none')
                #pylab.errorbar(rlog_r, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=1, marker = 'o', color='red', markersize=10, label='Red GC', linestyle='none')
                #pylab.errorbar(rlog_b, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=1,marker = 'o', color='blue', markersize=10, label='Blue GC', linestyle='none')        
                pylab.legend(loc='upper right', prop={'size':10})    
#                pylab.xlim(-1,1.5)
#                pylab.ylim(-4,2.5)
                pylab.xlabel('log(R)($arcmin$)')
                pylab.ylabel('log($\Sigma$)($N/arcmin^2$)')    
                pylab.show() 
        
        ################################################################################################################     
                
            ################################################################################################################    
            
            # decomment this part for adding IRAF/ELLIPSE tables data (please save in ASCII format 'ellipsetab.dat')            
            
            #rad = np.loadtxt('ellipsetab.dat', usecols=(1,))
            #rr = np.loadtxt('ellipsetab.dat', usecols=(2,))
            #err = np.loadtxt('ellipsetab.dat', usecols=(3,))
            #logrr = -2.5*np.log10(rr)
            #logerr = -2.5*0.4343*err/rr            
            
            plt = pylab.subplot(3,1,1)
            #pylab.plot(rad, logrr, marker='D', linestyle='none', markersize=5, color='red', label='Ellipse Fit', mfc='None')
            pylab.plot(x/60, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
            pylab.plot(x/60, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
            pylab.plot(x/60, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
            pylab.errorbar(rall/60, denslog+scale, yerr=errlog, capsize=6, capthick=0.5, ecolor='magenta', marker = 'o', markeredgecolor='magenta', markersize=10, mfc='None', label='All GC', linestyle='none')
            pylab.legend(loc='upper right', prop={'size':10})
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            plt2 = pylab.subplot(3,1,2)    
            pylab.plot(x/60, srsclog, color='red', linestyle='--')
            pylab.plot(x/60, exp1log, color ='blue', linestyle='--')
            pylab.plot(x/60, totallog, color='gray', linestyle='-')
            pylab.errorbar(rreds/60, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=0.5, ecolor='red', marker = 'o', markeredgecolor='red', markersize=10, mfc='None', label='Red GC', linestyle='none')
            pylab.legend(loc='upper right', prop={'size':10})    
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            plt3 = pylab.subplot(3,1,3)
            pylab.plot(x/60, srsclog, color='red', linestyle='--')
            pylab.plot(x/60, exp1log, color ='blue', linestyle='--')
            pylab.plot(x/60, totallog, color='gray', linestyle='-')
            pylab.errorbar(rblues/60, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=0.5, ecolor='blue',marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='Blue GC', linestyle='none')    
            pylab.legend(loc='upper right', prop={'size':10})    
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            pylab.subplots_adjust(hspace=.0001)
            pylab.setp(plt.get_xticklabels(), visible=False)
            pylab.setp(plt2.get_xticklabels(), visible=False)
            plt2.set_ylabel('log($\Sigma$) ($N/arcmin^2$)', fontsize=14)
            plt3.set_xlabel('(R) ($arcmin$)', fontsize=14)
            pylab.show()
        
    else:
        #### regular plots
        if c_sep == -1000:            #for the case where there is no bimodality
            N, area, rall = read_catalog_log(GCdensfile)     
            denslog, errlog, rlog = config_logplots(N, area, rall)    
            
            xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)    
            
            if op2 == 'y':
                PNedensfile = 'PNEdensity-circularbins.dat'
            
                dens, denserr, rmed, binsize = read_catalog(PNedensfile)                
                
                Npne, areapne, rpne = read_catalog_log(PNedensfile)
                denslogpne, errlogpne, rlogpne = config_logplots(Npne, areapne, rpne) 
                
                plt = pylab.subplot(1,1,1)
                pylab.plot(x, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
                pylab.plot(x, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
                pylab.plot(x, totallog, color='gray', linestyle='-', label = 'Light Profile from GALFIT')
                pylab.errorbar(rpne, denslogpne+scale, yerr=errlogpne, capsize=6, capthick=1, ecolor='green', marker = 'o', markeredgecolor='green', markersize=10, mfc='None', label='GMOS Sample', linestyle='none')
                pylab.errorbar(rall, denslog+scale, yerr=errlog, capsize=6, capthick=1, ecolor='blue', marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='WFC3 Sample', linestyle='none')                
                pylab.legend(loc='lower left', prop={'size':10})
                #pylab.xlim(-1,1)
                #pylab.ylim(-4,2.5)
                plt.set_ylabel('log($\Sigma$) ($N/arcsec^2$)', fontsize=14)
                plt.set_xlabel('R ($arcsec$)', fontsize=14)
                pylab.show()
            
            else:     
                
                
                srsc, exp1, total = galfit(x, mue, mud, bn, re, n, h)
                xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)
                
                #pylab.plot(rad, logrr, marker='D', linestyle='none', markersize=5, color='red', label='Ellipse Fit', mfc='None')
                plt = pylab.subplot(1,1,1)
                pylab.plot(x, srsc, color='red', label='Sersic Profile')
                pylab.plot(x, exp1, color ='blue', label='Exp Disk Profile')
                pylab.plot(x, total, color='gray', linestyle='--', label = 'Total Light Profile')
                pylab.errorbar(rall, denslog+scale, xerr=binsize_gc, yerr=GCdenserr, marker = 'o',color='purple', markersize=10, mfc='None', label='total GC density', linestyle='none')
                pylab.legend(loc='upper right', prop={'size':10})
                pylab.xlim(0,1100)
                #pylab.ylim(-10,10)
                pylab.xlabel('R($arcsec$)')
                pylab.ylabel('log($\Sigma$)($N/arcsec^2$)')                 
                pylab.show()
        else:
            if op2 == 'y':
        
                #optional: PNE
        #PNe 

                PNedensfile = 'PNEdensity-circularbins.dat'
            
                dens, denserr, rmed, binsize = read_catalog(PNedensfile)
                Npne, areapne, rpne = read_catalog_log(PNedensfile)
                denslogpne, errlogpne, rlogpne = config_logplots(Npne, areapne, rpne)
                srsc, exp1, total = galfit(x, mue, mud, bn, re, n, h)
                
                
                    
                    
                plt_all = pylab.subplot(1,1,1)
                pylab.plot(x, srsc, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
                pylab.plot(x, exp1, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
                pylab.plot(x, total, color='gray', linestyle='-', label = 'Total Light Profile')
                pylab.errorbar(rpne, denslogpne+scale, yerr=errlogpne, capsize=6, capthick=1, ecolor='green', marker = 'o', markeredgecolor='green', markersize=10, mfc='None', label='GMOS Sample', linestyle='none')
                #pylab.errorbar(rlog, denslog+scale, yerr=errlog, capsize=6, capthick=1, marker = 'o', color='magenta', markersize=10, label='All GC', linestyle='none')
                #pylab.errorbar(rlog_r, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=1, marker = 'o', color='red', markersize=10, label='Red GC', linestyle='none')
                #pylab.errorbar(rlog_b, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=1,marker = 'o', color='blue', markersize=10, label='Blue GC', linestyle='none')        
                pylab.legend(loc='lower left', prop={'size':10})    
                #pylab.xlim(-1,1.5)
                #pylab.ylim(-4,2.5)
                pylab.xlabel('R($arcsec$)')
                pylab.ylabel('log($\Sigma$)($N/arcsec^2$)')   
                #pylab.plt.gca().invert_yaxis()
                pylab.show() 
            
                
            N, area, rall = read_catalog_log(GCdensfile)
            N_r, area_r, rreds = read_catalog_log(redGCdensfile)
            N_b, area_b, rblues = read_catalog_log(blueGCdensfile)       
            denslog, errlog, rlog = config_logplots(N, area, rall)    
            denslog_r, errlog_r, rlog_r = config_logplots(N_r, area_r, rreds)
            denslog_b, errlog_b, rlog_b = config_logplots(N_b, area_b, rblues)
            
            xlog, srsclog, exp1log, totallog = log_galfit(x, mue, mud, bn, re, n, h)    
            
            plt = pylab.subplot(3,1,1)
            #pylab.plot(rad, logrr, marker='D', linestyle='none', markersize=5, color='red', label='Ellipse Fit', mfc='None')
            pylab.plot(x, srsclog, color='red', linestyle='--', label='Sersic Profile (Spheroid Light)')
            pylab.plot(x, exp1log, color ='blue', linestyle='--', label='Exp Disk Profile (Disk Light)')
            pylab.plot(x, totallog, color='gray', linestyle='-', label = 'Total Light Profile')
            pylab.errorbar(rall, denslog+scale, yerr=errlog, capsize=6, capthick=0.5, ecolor='magenta', marker = 'o', markeredgecolor='magenta', markersize=10, mfc='None', label='All GC', linestyle='none')
            pylab.legend(loc='upper right', prop={'size':10})
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            plt2 = pylab.subplot(3,1,2)    
            pylab.plot(x, srsclog, color='red', linestyle='--')
            pylab.plot(x, exp1log, color ='blue', linestyle='--')
            pylab.plot(x, totallog, color='gray', linestyle='-')
            pylab.errorbar(rreds, denslog_r+scale, yerr=errlog_r, capsize=6, capthick=0.5, ecolor='red', marker = 'o', markeredgecolor='red', markersize=10, mfc='None', label='Red GC', linestyle='none')
            pylab.legend(loc='upper right', prop={'size':10})    
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            plt3 = pylab.subplot(3,1,3)
            pylab.plot(x, srsclog, color='red', linestyle='--')
            pylab.plot(x, exp1log, color ='blue', linestyle='--')
            pylab.plot(x, totallog, color='gray', linestyle='-')
            pylab.errorbar(rblues, denslog_b+scale, yerr=errlog_b, capsize=6, capthick=0.5, ecolor='blue',marker = 'o', markeredgecolor='blue', markersize=10, mfc='None', label='Blue GC', linestyle='none')    
            pylab.legend(loc='upper right', prop={'size':10})    
            #pylab.xlim(-0.23,1.38)
            #pylab.ylim(-4,2.5)
            
            pylab.subplots_adjust(hspace=.0001)
            pylab.setp(plt.get_xticklabels(), visible=False)
            pylab.setp(plt2.get_xticklabels(), visible=False)
            plt2.set_ylabel('log($\Sigma$) ($N/arcmin^2$)', fontsize=14)
            plt3.set_xlabel('(R) ($arcmin$)', fontsize=14)
            pylab.show()
    