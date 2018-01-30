import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

with open('current_galaxy.dat', 'r') as f:     
    gal = [line.split()[0] for line in f]
    gal = gal[0] 
    
op = raw_input('is this galaxy unimodal? y/n ')

def load_like_uni2(op):
        
    V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
    Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
    Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
    sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
    sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
    sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
    R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
    Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
    Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
    Rmin = R - Rmin
    Rmax = Rmax - R
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax            

def load_like_bi2(op):
    
    if op == 'allGC':
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R   
        
    if op=='redGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
         
        
    if op=='blueGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
                
    #binsizes
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax

def load_like_uni(op):
    if op=='GC':        
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R
                
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V
    
        sigmaerr_min = sigma - sigmaerr_min
        sigmaerr_max = sigmaerr_max - sigma
        
    if op=='PNe':
        V = np.loadtxt('pnelike_'+gal+'.dat', usecols=(1,))
        Verr_min = np.loadtxt('pneerr_'+gal+'.dat', usecols=(0,))
        Verr_max = np.loadtxt('pneerr_'+gal+'.dat', usecols=(1,))
        sigma = np.loadtxt('pnelike_'+gal+'.dat', usecols=(2,))
        sigmaerr_min = np.loadtxt('pneerr_'+gal+'.dat', usecols=(2,))
        sigmaerr_max = np.loadtxt('pneerr_'+gal+'.dat', usecols=(3,))
        R = np.loadtxt('pnelike_'+gal+'.dat', usecols=(7,))
        
        Rmin = np.loadtxt('pnebin_'+gal+'.dat', usecols=(3,))
        Rmax = np.loadtxt('pnebin_'+gal+'.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V
    
        sigmaerr_min = sigma - sigmaerr_min
        sigmaerr_max = sigmaerr_max - sigma
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax
    
def load_like_bi(op):
    
    if op == 'allGC':
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V

        sigmaerr_min = sigma - sigmaerr_min
        sigmaerr_max = sigmaerr_max - sigma    
        
    if op=='redGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_red.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_red.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_red.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V        

        sigmaerr_min = sigma - sigmaerr_min
        sigmaerr_max = sigmaerr_max - sigma  
        
    if op=='blueGC':    
        V = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(10,))
        Verr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(10,))
        Verr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(11,))
        sigma = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(4,))
        sigmaerr_min = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(6,))
        sigmaerr_max = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/error_blue.dat', usecols=(7,))
        R = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/like_blue.dat', usecols=(7,))

        Rmin = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/2.5/N'+gal+'/bin_blue.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R  
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V        
        
        sigmaerr_min = sigma - sigmaerr_min
        sigmaerr_max = sigmaerr_max - sigma          
    #binsizes
    
    return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax

def load_twocompPNE():
        V = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(1,))
        Verr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(0,))
        Verr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(1,))
        sigma_r = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(2,))
        sigmarerr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(2,))
        sigmarerr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(3,))
        sigma_p = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(3,))
        sigmaperr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(4,))
        sigmaperr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(5,))
        R = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(7,))
        
        Vb = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(10,))
        Vb_err_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(10,))
        Vb_err_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(11,))     
        sigma_b = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(4,))
        sigma_b_err_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(6,))
        sigma_b_err_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/pne_likes/bin_'+gal+'.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/pne_likes/bin_'+gal+'.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V
        Vb_err_min = Vb - Vb_err_min
        Vb_err_max = Vb_err_max - Vb
        sigma_b_err_min = sigma_b - sigma_b_err_min
        sigma_b_err_max = sigma_b_err_max - sigma_b
        sigmarerr_min = sigma_r - sigmarerr_min
        sigmarerr_max = sigmarerr_max - sigma_r
        sigmaperr_min = sigma_p - sigmaperr_min
        sigmaperr_max = sigmaperr_max - sigma_p

        sigma = (sigma_r+sigma_p)/2
        sigmaerr_max = (sigmarerr_max+sigmaperr_max)/2
        sigmaerr_min = (sigmarerr_min+sigmaperr_min)/2
        
        return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax, Vb, Vb_err_min, Vb_err_max, sigma_b, sigma_b_err_min, sigma_b_err_max
        
def load_twocompPNE2():
        V = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(1,))
        Verr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(0,))
        Verr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(1,))
        sigma_r = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(2,))
        sigmarerr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(2,))
        sigmarerr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(3,))
        sigma_p = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(3,))
        sigmaperr_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(4,))
        sigmaperr_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(5,))
        R = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(7,))
            
        sigma_b = np.loadtxt('/home/emilio/MLE/pne_likes/likelihood_'+gal+'.dat', usecols=(4,))
        sigma_b_err_min = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(6,))
        sigma_b_err_max = np.loadtxt('/home/emilio/MLE/pne_likes/error_'+gal+'.dat', usecols=(7,))
        
        Rmin = np.loadtxt('/home/emilio/MLE/pne_likes/bin_'+gal+'.dat', usecols=(3,))
        Rmax = np.loadtxt('/home/emilio/MLE/pne_likes/bin_'+gal+'.dat', usecols=(4,))
        Rmin = R - Rmin
        Rmax = Rmax - R
        
        Verr_min = V - Verr_min
        Verr_max = Verr_max - V
        sigma_b_err_min = sigma_b - sigma_b_err_min
        sigma_b_err_max = sigma_b_err_max - sigma_b
        sigmarerr_min = sigma_r - sigmarerr_min
        sigmarerr_max = sigmarerr_max - sigma_r
        sigmaperr_min = sigma_p - sigmaperr_min
        sigmaperr_max = sigmaperr_max - sigma_p

        sigma = (sigma_r+sigma_p)/2
        sigmaerr_max = (sigmarerr_max+sigmaperr_max)/2
        sigmaerr_min = (sigmarerr_min+sigmaperr_min)/2
        
        return V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax, sigma_b, sigma_b_err_min, sigma_b_err_max

if op=='y':
#    V = np.loadtxt('like_.dat', usecols=(10,))
#    Verr_min = np.loadtxt('error_.dat', usecols=(10,))
#    Verr_max = np.loadtxt('error_.dat', usecols=(11,))
#    sigma = np.loadtxt('like_.dat', usecols=(4,))
#    sigmaerr_min = np.loadtxt('error_.dat', usecols=(6,))
#    sigmaerr_max = np.loadtxt('error_.dat', usecols=(7,))
#    R = np.loadtxt('like_.dat', usecols=(7,))
#    
#    Rmin = np.loadtxt('bin_.dat', usecols=(3,))
#    Rmax = np.loadtxt('bin_.dat', usecols=(4,))
#    Rmin = R - Rmin
#    Rmax = Rmax - R
#    
#    Verr_min = V - Verr_min
#    Verr_max = Verr_max - V
#
#    sigmaerr_min = sigma - sigmaerr_min
#    sigmaerr_max = sigmaerr_max - sigma

    gV, gVerr_min, gVerr_max, gsigma, gsigmaerr_min, gsigmaerr_max, gR, gRmin, gRmax = load_like_uni('GC')
    pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax = load_like_uni('PNe')   
    
    plot1 = plt.subplot(2, 1, 1)
    plt.errorbar(pR, pV, marker='o', xerr=[pRmin, pRmax], yerr=[pVerr_min, pVerr_max], color='green', linestyle='None', label='PNe')    
    plt.errorbar(gR, gV, marker='o', xerr=[gRmin, gRmax], yerr=[gVerr_min, gVerr_max], color='purple', linestyle='None', label='All GCs')
    plt.ylabel('V (km/s)')
    plt.ylim(0, 200)
    plt.xlim(0, 150)
    plt.legend(loc='upper left')
    
    plot4 = plt.subplot(2, 1, 2)
    plt.errorbar(pR, psigma ,marker='o', xerr=[pRmin, pRmax], yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None')
    plt.errorbar(gR, gsigma ,marker='o', xerr=[gRmin, gRmax], yerr=[gsigmaerr_min, gsigmaerr_max], color='purple', linestyle='None')
    plt.ylabel('$\sigma$ (km/s)')
    plt.xlabel('R ($arcsec$)')
    plt.ylim(0, 90)
    plt.xlim(0, 150)    
    
    loc = plticker.MultipleLocator(base = 50.0)
    loc2 = plticker.MultipleLocator(base = 50.0)
    loc3 = plticker.MultipleLocator(base = 20.0)    
    plot1.yaxis.set_major_locator(loc)
    plot4.xaxis.set_major_locator(loc)
    plot1.xaxis.set_major_locator(loc2)
    plot4.yaxis.set_major_locator(loc3)
    plt.setp(plot1.get_xticklabels(), visible=False)
    
    plt.subplots_adjust(hspace=.0001)
    plt.show()
    
    
######################################################33 new plot

    pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax = load_like_uni('PNe')
       
    
    p1 = plt.subplot(2,1,1)
        
    plt.errorbar(pR, pV, marker='o', yerr=[pVerr_min, pVerr_max], color='green', linestyle='None', markersize=10)
    plt.errorbar(gR, gV,marker='o', yerr=[gVerr_min, gVerr_max], markeredgewidth=1,markeredgecolor='purple', markerfacecolor='None', linestyle='None', label='All GCs', markersize=12, ecolor='purple')
    plt.ylabel('$V_{Disk}$')
    plt.ylim(0, 300)
    plt.xlim(30, 275)


    p2 = plt.subplot(2,1,2)

    plt.errorbar(pR, psigma, marker='o', yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None', markersize=10)
    plt.errorbar(gR, gsigma,marker='o', yerr=[gsigmaerr_min, gsigmaerr_max], markeredgewidth=1,markeredgecolor='purple', markerfacecolor='None', linestyle='None', label='All GCs', markersize=12, ecolor='purple')
    plt.ylabel('$\sigma_{disk}$')
    plt.ylim(0, 300)
    plt.xlim(30, 275)
    
     
    p1.yaxis.set_major_locator(loc3)    
    p2.yaxis.set_major_locator(loc3)    
    p2.set_yticks((0, 100, 200), minor=False)	
    plt.setp(p1.get_xticklabels(), visible=False)    
    plt.subplots_adjust(hspace=.0001)
    plt.show()   
    
    
    #    TODO: Correct this errorbars here. They are very, very weird...        
    V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax = load_like_uni2('allGC')
#    pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax, pVb, pVb_err_min, pVb_err_max, psigma_b, psigma_b_err_min, psigma_b_err_max = load_twocompPNE()
   
    Vos = V/sigma

    Vos_max = Verr_max/sigmaerr_max
    Vos_min = Verr_min/sigmaerr_min
    print 'V/sigma for all: ', np.mean(Vos), np.mean(Vos_max), np.mean(Vos_min)
    #how to proceed with error propagation with assymetrical errorbars? let's proceed without it for now:
    
    plt.errorbar(R, Vos, yerr=[Vos_min, Vos_max], markerfacecolor='None', markersize=10, ecolor='purple', markeredgecolor='purple', marker='o', linestyle='none')
    plt.show()
else:
#    V = np.loadtxt('like_.dat', usecols=(10,))
#    Verr_min = np.loadtxt('error_.dat', usecols=(10,))
#    Verr_max = np.loadtxt('error_.dat', usecols=(11,))
#    sigma = np.loadtxt('like_.dat', usecols=(4,))
#    sigmaerr_min = np.loadtxt('error_.dat', usecols=(6,))
#    sigmaerr_max = np.loadtxt('error_.dat', usecols=(7,))
#    R = np.loadtxt('like_.dat', usecols=(7,))
#    
#    V_r = np.loadtxt('like_red.dat', usecols=(10,))
#    Verr_r_min = np.loadtxt('error_red.dat', usecols=(10,))
#    Verr_r_max = np.loadtxt('error_red.dat', usecols=(11,))
#    sigma_r = np.loadtxt('like_red.dat', usecols=(4,))
#    sigmaerr_r_min = np.loadtxt('error_red.dat', usecols=(6,))
#    sigmaerr_r_max = np.loadtxt('error_red.dat', usecols=(7,))
#    R_r = np.loadtxt('like_red.dat', usecols=(7,))
#    
#    V_b = np.loadtxt('like_blue.dat', usecols=(10,))
#    Verr_b_min = np.loadtxt('error_blue.dat', usecols=(10,))
#    Verr_b_max = np.loadtxt('error_blue.dat', usecols=(11,))
#    sigma_b = np.loadtxt('like_blue.dat', usecols=(4,))
#    sigmaerr_b_min = np.loadtxt('error_blue.dat', usecols=(6,))
#    sigmaerr_b_max = np.loadtxt('error_blue.dat', usecols=(7,))
#    R_b = np.loadtxt('like_blue.dat', usecols=(7,))
#    
#    #binsizes
#    
#    Rmin = np.loadtxt('bin_.dat', usecols=(3,))
#    Rmax = np.loadtxt('bin_.dat', usecols=(4,))
#    Rmin = R - Rmin
#    Rmax = Rmax - R
#    
#    Rmin_r = np.loadtxt('bin_red.dat', usecols=(3,))
#    Rmax_r = np.loadtxt('bin_red.dat', usecols=(4,))
#    Rmin_r = R_r - Rmin_r
#    Rmax_r = Rmax_r - R_r
#    
#    Rmin_b = np.loadtxt('bin_blue.dat', usecols=(3,))
#    Rmax_b = np.loadtxt('bin_blue.dat', usecols=(4,))
#    Rmin_b = R_b - Rmin_b
#    Rmax_b = Rmax_b - R_b
#    
#    
#    Verr_min = V - Verr_min
#    Verr_max = Verr_max - V
#    Verr_r_min = V_r - Verr_r_min
#    Verr_r_max = Verr_r_max - V_r
#    Verr_b_min = V_b - Verr_b_min
#    Verr_b_max = Verr_b_max - V_b
#    
#    sigmaerr_min = sigma - sigmaerr_min
#    sigmaerr_max = sigmaerr_max - sigma
#    sigmaerr_r_min = sigma_r - sigmaerr_r_min
#    sigmaerr_r_max = sigmaerr_r_max - sigma_r
#    sigmaerr_b_min = sigma_b - sigmaerr_b_min
#    sigmaerr_b_max = sigmaerr_b_max - sigma_b
    
    #convert all things to arcmin? 
    
    #R = R/60
    #R_r = R_r/60
    #R_b = R_b/60
    #Rmin = Rmin/60
    #Rmax = Rmax/60
    #Rmin_r = Rmin_r/60
    #Rmax_r = Rmax_r/60
    #Rmin_b = Rmin_b/60
    #Rmax_b = Rmax_b/60
    
    V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax = load_like_bi('allGC')
    V_r, Verr_r_min, Verr_r_max, sigma_r, sigmaerr_r_min, sigmaerr_r_max, R_r, Rmin_r, Rmax_r = load_like_bi('redGC')
    V_b, Verr_b_min, Verr_b_max, sigma_b, sigmaerr_b_min, sigmaerr_b_max, R_b, Rmin_b, Rmax_b = load_like_bi('blueGC')
    pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax = load_like_uni('PNe')
    
    plot1 = plt.subplot(2, 3, 1)
    plt.errorbar(pR, pV, marker='o', xerr=[pRmin, pRmax], yerr=[pVerr_min, pVerr_max], color='green', linestyle='None', label='PNe') 
    plt.errorbar(R, V, marker='o', xerr=[Rmin, Rmax], yerr=[Verr_min, Verr_max], color='purple', linestyle='None', label='All GCs')
    plt.ylabel('V (km/s)')
    plt.ylim(0, 450)
    plt.legend(loc='upper left')
    plot2 = plt.subplot(2, 3, 2)
    #plt.title('NGC 2768')
    plt.errorbar(pR, pV, marker='o', xerr=[pRmin, pRmax], yerr=[pVerr_min, pVerr_max], color='green', linestyle='None') 
    plt.errorbar(R_r, V_r,marker='o', xerr=[Rmin_r, Rmax_r], yerr=[Verr_r_min, Verr_r_max], color='red', linestyle='None', label='Red GCs')
    plt.ylim(0, 450)
    plt.legend(loc='upper left')
    plot3 = plt.subplot(2, 3, 3)
    plt.errorbar(pR, pV, marker='o', xerr=[pRmin, pRmax], yerr=[pVerr_min, pVerr_max], color='green', linestyle='None') 
    plt.errorbar(R_b, V_b,marker='o', xerr=[Rmin_b, Rmax_b], yerr=[Verr_b_min, Verr_b_max], color='blue', linestyle='None', label='Blue GCs')
    plt.ylim(0, 450)
    plt.legend(loc='upper left')
    plot4 = plt.subplot(2, 3, 4)
    plt.errorbar(pR, psigma ,marker='o', xerr=[pRmin, pRmax], yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None')
    plt.errorbar(R, sigma,marker='o', xerr=[Rmin, Rmax], yerr=[sigmaerr_min, sigmaerr_max], color='purple', linestyle='None')
    plt.ylabel('$\sigma$ (km/s)')
    plt.ylim(0, 450)
    plot5 = plt.subplot(2, 3, 5)
    plt.errorbar(pR, psigma ,marker='o', xerr=[pRmin, pRmax], yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None')
    plt.errorbar(R_r, sigma_r,marker='o', xerr=[Rmin_r, Rmax_r], yerr=[sigmaerr_r_min, sigmaerr_r_max], color='red', linestyle='None')
    plt.xlabel('R (arcsec)')
    plt.ylim(0, 450)
    plot6 = plt.subplot(2, 3, 6)
    plt.errorbar(pR, psigma ,marker='o', xerr=[pRmin, pRmax], yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None')
    plt.errorbar(R_b, sigma_b,marker='o', xerr=[Rmin_b, Rmax_b], yerr=[sigmaerr_b_min, sigmaerr_b_max], color='blue', linestyle='None')
    plt.ylim(0, 450)
    
    loc = plticker.MultipleLocator(base = 200.0)
    loc2 = plticker.MultipleLocator(base=100.0)
    loc3 = plticker.MultipleLocator(base = 100.0)
    plot1.yaxis.set_major_locator(loc3)
    plot4.xaxis.set_major_locator(loc)
    plot1.xaxis.set_major_locator(loc)
    plot4.yaxis.set_major_locator(loc3)
    plot5.yaxis.set_major_locator(loc3)
    plot2.yaxis.set_major_locator(loc3)
    plot3.yaxis.set_major_locator(loc3)
    plot6.yaxis.set_major_locator(loc3)
    plot5.xaxis.set_major_locator(loc)
    plot2.xaxis.set_major_locator(loc)
    plot6.xaxis.set_major_locator(loc)
    plot3.xaxis.set_major_locator(loc)
    
    plt.setp(plot5.get_yticklabels(), visible=False)
    plt.setp(plot6.get_yticklabels(), visible=False)
    plt.setp(plot1.get_xticklabels(), visible=False)
    plt.setp(plot2.get_xticklabels(), visible=False)
    plt.setp(plot3.get_xticklabels(), visible=False)
    plt.setp(plot2.get_yticklabels(), visible=False)
    plt.setp(plot3.get_yticklabels(), visible=False)
    
    plt.subplots_adjust(hspace=.0001)
    plt.subplots_adjust(wspace=.0001)
    plt.show()
########################################## new plot
    
    if gal == '2768':
    
        pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax, pVb, pVb_err_min, pVb_err_max, psigma_b, psigma_b_err_min, psigma_b_err_max = load_twocompPNE()
        
        p1 = plt.subplot(4,1,1)
        
        plt.errorbar(pR, pV, marker='o', yerr=[pVerr_min, pVerr_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, V_r,marker='o', yerr=[Verr_r_min, Verr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, V_b,marker='o', yerr=[Verr_b_min, Verr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$V_{Disk}$')
        plt.ylim(0, 300)
        plt.xlim(30, 275)
    
        
        p2 = plt.subplot(4,1,2)
        
        plt.errorbar(pR, psigma, marker='o', yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, sigma_r,marker='o', yerr=[sigmaerr_r_min, sigmaerr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, sigma_b,marker='o', yerr=[sigmaerr_b_min, sigmaerr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$\sigma_{disk}$')
        plt.ylim(0, 300)
        plt.xlim(30, 275)        
                
        p3 = plt.subplot(4,1,3)
        
        plt.errorbar(pR, pVb, marker='o', yerr=[pVb_err_min, pVb_err_max], color='green', linestyle='None', label='PNe', markersize=10, ecolor='g')
        plt.errorbar(R_r, V_r,marker='o', yerr=[Verr_r_min, Verr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, V_b,marker='o', yerr=[Verr_b_min, Verr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$V_{Sph}$')
        plt.ylim(0, 300)
        plt.xlim(30, 275)
        
        p4 = plt.subplot(4,1,4)
        
        plt.errorbar(pR, psigma_b, marker='o', yerr=[psigma_b_err_min, psigma_b_err_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, sigma_r,marker='o', yerr=[sigmaerr_r_min, sigmaerr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, sigma_b,marker='o', yerr=[sigmaerr_b_min, sigmaerr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$\sigma_{Sph}$')
        plt.ylim(0, 300)
        plt.xlim(30, 275)
    
     
        p1.yaxis.set_major_locator(loc3)    
        p2.yaxis.set_major_locator(loc3)    
        p3.yaxis.set_major_locator(loc3)
        p4.yaxis.set_major_locator(loc3)
        p2.set_yticks((0, 100, 200), minor=False)	
        p3.set_yticks((0, 100, 200), minor=False)	
        p4.set_yticks((0, 100, 200), minor=False)	
        plt.setp(p1.get_xticklabels(), visible=False)    
        plt.setp(p2.get_xticklabels(), visible=False)    
        plt.setp(p3.get_xticklabels(), visible=False)
        plt.subplots_adjust(hspace=.0001)
        plt.show()
    
    else:
        
        #TODO: Add sigma_sph to this galaxies and use the average of sigmas for sigma_disk (I think I already done this)    
        
        pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax, psigma_b, psigma_b_err_min, psigma_b_err_max = load_twocompPNE2()       
       
    
        p1 = plt.subplot(3,1,1)
        
        plt.errorbar(pR, pV, marker='o', yerr=[pVerr_min, pVerr_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, V_r,marker='o', yerr=[Verr_r_min, Verr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, V_b,marker='o', yerr=[Verr_b_min, Verr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$V_{Disk}$')
        plt.ylim(0, 300)
        plt.xlim(7, 381)
    
        
        p2 = plt.subplot(3,1,2)
        
        plt.errorbar(pR, psigma, marker='o', yerr=[psigmaerr_min, psigmaerr_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, sigma_r,marker='o', yerr=[sigmaerr_r_min, sigmaerr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, sigma_b,marker='o', yerr=[sigmaerr_b_min, sigmaerr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$\sigma_{disk}$')
        plt.ylim(0, 300)
        plt.xlim(7, 381)
        
        p3 = plt.subplot(3,1,3)
        
        plt.errorbar(pR, psigma_b, marker='o', yerr=[psigma_b_err_min, psigma_b_err_max], color='green', linestyle='None', markersize=10)
        plt.errorbar(R_r, sigma_r,marker='o', yerr=[sigmaerr_r_min, sigmaerr_r_max], markeredgewidth=1,markeredgecolor='red', markerfacecolor='None', linestyle='None', label='Red GCs', markersize=12, ecolor='red')
        plt.errorbar(R_b, sigma_b,marker='o', yerr=[sigmaerr_b_min, sigmaerr_b_max], markeredgewidth=1,markeredgecolor='blue', markerfacecolor='None', linestyle='None', label='Blue GCs', markersize=12, ecolor='blue')
        plt.ylabel('$\sigma_{Sph}$')
        plt.ylim(0, 300)
        plt.xlim(7, 381)
    
     
        p1.yaxis.set_major_locator(loc3)    
        p2.yaxis.set_major_locator(loc3)  
        p3.yaxis.set_major_locator(loc3)
        p2.set_yticks((0, 100, 200), minor=False)	
        p3.set_yticks((0, 100, 200), minor=False)
        plt.setp(p1.get_xticklabels(), visible=False)    
        plt.setp(p2.get_xticklabels(), visible=False)  
        plt.subplots_adjust(hspace=.0001)
        plt.show()        
        
################ another new plot
#    TODO: Correct this errorbars here. They are very, very weird...        
    V, Verr_min, Verr_max, sigma, sigmaerr_min, sigmaerr_max, R, Rmin, Rmax = load_like_bi2('allGC')
    V_r, Verr_r_min, Verr_r_max, sigma_r, sigmaerr_r_min, sigmaerr_r_max, R_r, Rmin_r, Rmax_r = load_like_bi2('redGC')
    V_b, Verr_b_min, Verr_b_max, sigma_b, sigmaerr_b_min, sigmaerr_b_max, R_b, Rmin_b, Rmax_b = load_like_bi2('blueGC')
#    pV, pVerr_min, pVerr_max, psigma, psigmaerr_min, psigmaerr_max, pR, pRmin, pRmax, pVb, pVb_err_min, pVb_err_max, psigma_b, psigma_b_err_min, psigma_b_err_max = load_twocompPNE()
   
    Vos = V/sigma
    Vosreds = V_r/sigma_r
    Vosblues = V_b/sigma_b
    
    Vos_max = Verr_max/sigmaerr_max
    Vos_min = Verr_min/sigmaerr_min
    Vosreds_max = Verr_r_max/sigmaerr_r_max
    Vosreds_min = Verr_r_min/sigmaerr_r_min
    Vosblues_max = Verr_b_max/sigmaerr_b_max
    Vosblues_min = Verr_b_min/sigmaerr_b_min
    
    print 'V/sigma for all: ', np.mean(Vos), np.mean(Vos_max), np.mean(Vos_min)
    print 'V/sigma for reds: ', np.mean(Vosreds), np.mean(Vosreds_max), np.mean(Vosreds_min)
    print 'V/sigma for blues: ', np.mean(Vosblues), np.mean(Vosblues_max), np.mean(Vosblues_min)
    
    #how to proceed with error propagation with assymetrical errorbars? let's proceed without it for now:
    
    plt.errorbar(R, Vos, yerr=[Vos_min, Vos_max], markerfacecolor='None', markersize=10, ecolor='purple', markeredgecolor='purple', marker='o', linestyle='none')
    plt.errorbar(R_r, Vosreds, yerr=[Vosreds_min, Vosreds_max], markerfacecolor='None',ecolor='red', markersize=10,markeredgecolor='red', marker='o', linestyle='none')
    plt.errorbar(R_b, Vosblues, yerr=[Vosblues_min, Vosblues_max], markerfacecolor='None',ecolor='blue', markersize=10,markeredgecolor='blue', marker='o', linestyle='none')
    plt.show()
    
    
    
