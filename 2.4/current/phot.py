from pyraf import iraf
import os
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.pyplot import xticks
import ast

#############################################################################################
##Pyraf script to create the fmap. Inputs: galfit bulge and total model fits, RA&DEC catalog.
##Outputs: inputs to prob_bulge_bd.f and fmap plot.  
##TO DO: > add some way to run galfit directly, given a galfit input file. (27/10/17: is this really necessary?)
##       > fix the problems with coordinates of galcenter for some galaxies
#############################################################################################

def gen_cat(RA, DEC, gal_center_x, gal_center_y):
     DECdeg = DEC*u.degree
     DECrad = DECdeg.to(u.rad)
     DECrad = DECrad.value
     for i in range(0, len(RA)):
		RA[i] = RA[i]*3600
		DEC[i] = DEC[i]*3600  
     path = '/home/emilio/MLE/Galaxies/'+gal+'/'  
     galinput = 'N'+gal+'input.dat'
     galinput = path+galinput     
     inp = [x.split(' ')[0] for x in open(galinput).readlines()]     
     RAgal = inp[0]
     DECgal = inp[1]                     
     galcenter = SkyCoord(RAgal, DECgal)     #real galaxy center
     RAg = galcenter.ra.arcsec
     DECg = galcenter.dec.arcsec
     
     print('Galaxy center x: ', gal_center_x)
     print('Galaxy center y: ', gal_center_y)     
     
     
     DRA = []
     DDEC = []
     
     galcentery = (galcenter.dec.degree)
     galcentery = np.pi*galcentery/180
     
     print galcenter.dec.degree
     print galcentery
     for i in range(0, len(RA)):
		DRA.append((-RA[i]+RAg)*np.cos(galcentery))
		DDEC.append((DEC[i]-DECg))
     for i in range(0, len(RA)):
		DRA[i] = DRA[i] + gal_center_x
		DDEC[i] = DDEC[i] + gal_center_y
     catname = path+'2MASS.GC.cat'
     with open(catname, 'w') as fcat:
		for i in range(0, len(RA)):
			print>>fcat, int(DRA[i]), int(DDEC[i])

def rename_columns_catalogs(list_of_catalogs):                                   ######### Made by Loic Le Tiran
    """ Gives names to the catalog columns outputs by IRAF
    :param list_of_catalogs: a list of tables to rename columns in
    :return: the upgraded list
    """
    for t in list_of_catalogs:
        t.rename_column('col1', 'X')
        t.rename_column('col2', 'Y')
        t.rename_column('col3', 'error_position')
        t.rename_column('col4', 'flux')
        t.rename_column('col5', 'sum')
        t.rename_column('col6', 'area')
        t.rename_column('col7', 'ID')
    return list_of_catalogs

#def ffinder(catalog_bulge, catalog_total, catalog_fraction, cat_flux_per_area):   ########## Made by Loic Le Tiran
#    """ Gives fluxes of the PN per area and saves them in a file
#    :param catalog_bulge: the IRAF catalog with the flux values from the bulge
#    :param catalog_total: the IRAF catalog with the flux values from the total
#    :param catalog_fraction: the IRAF catalog with the flux values from the division
#    :param cat_flux_per_area: the output catalog with the fluxes per area
#    :return:
#    """
#    bulge = Table.read(catalog_bulge, format="ascii")
#    total = Table.read(catalog_total, format="ascii")
#    fraction = Table.read(catalog_fraction, format="ascii")
#    list_tables = [bulge, total, fraction]
#    list_tables = rename_columns_catalogs(list_tables)
#    # Replaces the nan values by 0.5
#    f_bulge_over_total = bulge["flux"] / total["flux"]
#    f_fraction_per_area = fraction["flux"] / fraction["area"]
#    f_bulge_over_total[np.isnan(f_bulge_over_total)] = .5
#    f_fraction_per_area[np.isnan(f_fraction_per_area)] = .5
#    f = open(cat_flux_per_area, "w")
#    for i in np.arange(len(total)):
#        f.write(str(f_bulge_over_total[i]) +" ")
#        f.write(str(f_fraction_per_area[i]) +" ")
#        f.write(str(bulge["ID"][i]) +"\n")
#    f.close()

def ffinder(catalog_fraction, cat_flux_per_area):   ########## Made by Loic Le Tiran
    """ Gives fluxes of the PN per area and saves them in a file
    :param catalog_bulge: the IRAF catalog with the flux values from the bulge
    :param catalog_total: the IRAF catalog with the flux values from the total
    :param catalog_fraction: the IRAF catalog with the flux values from the division
    :param cat_flux_per_area: the output catalog with the fluxes per area
    :return:
    """
    #bulge = Table.read(catalog_bulge, format="ascii")
    #total = Table.read(catalog_total, format="ascii")
    fraction = Table.read(catalog_fraction, format="ascii")
    list_tables = [fraction]
    list_tables = rename_columns_catalogs(list_tables)
    # Replaces the nan values by 0.5
    #f_bulge_over_total = bulge["flux"] / total["flux"]
    f_fraction_per_area = fraction["flux"] / fraction["area"]
    #f_bulge_over_total[np.isnan(f_bulge_over_total)] = .5
    f_fraction_per_area[np.isnan(f_fraction_per_area)] = .5
    f = open(cat_flux_per_area, "w")
    for i in np.arange(len(fraction)):
        #f.write(str(f_bulge_over_total[i]) +" ")
        print >>f, f_fraction_per_area[i]
        #f.write(str(bulge["ID"][i]) +"\n")
    f.close()
      

def plot_fmap(gal_center_x, gal_center_y):
     DRA = np.loadtxt(path+'2MASS.GC.cat', usecols=(0,))
     DDEC = np.loadtxt(path+'2MASS.GC.cat', usecols=(1,))
     f = np.loadtxt(path+'N'+gal+'_f.GC.dat', usecols=(0,))


     fits_file = path+'f.fits'
     hdu = fits.open(fits_file)[0]
     
     fig, ax = plt.subplots()
     plt.imshow(hdu.data, origin='lower',cmap='gray')
     plt.scatter(DRA, DDEC, c=f, cmap='Blues') 
     plt.xlim(102.16, 682.36)
     plt.ylim(158.99, 610.43)
     xlabels = ax.get_xticks().tolist()
     for i in range(0, len(xlabels)):
         xlabels[i] = np.around(xlabels[i]-gal_center_x, decimals=-1)
         if xlabels[i] == -0.0:
             xlabels[i] = 0.0   
     ax.set_xticklabels(xlabels)    
     ylabels = ax.get_yticks().tolist()
     for i in range(0, len(ylabels)):
         ylabels[i] = np.around(ylabels[i]-gal_center_y, decimals=-1)
         if ylabels[i] == -0.0:
             ylabels[i] = 0.0   
     ax.set_yticklabels(ylabels)
     ax.set_xlabel('$\Delta RA (arcsec)$', fontsize=20)
     ax.set_ylabel('$\Delta DEC (arcsec)$', fontsize=20)
     cb = plt.colorbar(fraction=0.04)
     cb.set_label('$f_{i}$', fontsize=20)
    	#plt.set_xlabel('RA')   
    	#plt.set_ylabel('DEC')     
     plt.show()
     
 
def plot_vmap(gal_center_x, gal_center_y):
    DRA = np.loadtxt(path+'2MASS.GC.cat', usecols=(0,))
    DDEC = np.loadtxt(path+'2MASS.GC.cat', usecols=(1,))
    v = np.loadtxt(path+'N'+gal+'GC.dat', usecols=(3,))


    fits_file = path+'f.fits'
    hdu = fits.open(fits_file)[0]

    Vmin = min(v)
    Vmax = max(v)
    
    fig, ax = plt.subplots()
    plt.imshow(hdu.data, origin='lower',cmap='gray')
    plt.scatter(DRA, DDEC, c=v, vmin=Vmin, vmax=Vmax)
#    plt.plot(gal_center_x, gal_center_y, 'ro')
    plt.xlim(102.16, 682.36)
    plt.ylim(158.99, 610.43)
    xlabels = ax.get_xticks().tolist()
    for i in range(0, len(xlabels)):
        xlabels[i] = np.around(xlabels[i]-gal_center_x, decimals=-1)
        if xlabels[i] == -0.0:
            xlabels[i] = 0.0   
    ax.set_xticklabels(xlabels)    
    ylabels = ax.get_yticks().tolist()
    for i in range(0, len(ylabels)):
        ylabels[i] = np.around(ylabels[i]-gal_center_y, decimals=-1)
        if ylabels[i] == -0.0:
            ylabels[i] = 0.0   
    ax.set_yticklabels(ylabels) 
    cb = plt.colorbar(fraction=0.04)
    cb.set_label('Radial Velocity (km/s)')
    plt.xlabel('')
    plt.show() 

def read_catalog_gc(galcat):
    RA = np.loadtxt(galcat, usecols=(1,))
    DEC = np.loadtxt(galcat, usecols=(2,))       
    V = np.loadtxt(galcat, usecols=(3,))            
    col1 = np.loadtxt(galcat, usecols=(4,))
    col2 = np.loadtxt(galcat, usecols=(5,))
    colour = col1 - col2    
    return RA, DEC, V, colour
    
def read_catalog_pne(pnecat):
    RAp = np.loadtxt(pnecat, usecols=(1,))
    DECp = np.loadtxt(pnecat, usecols=(2,))
    
    Vp = np.loadtxt(pnecat, usecols=(3,))
    
    return RAp, DECp, Vp
    
def plot_histf(ffile):
    import matplotlib.pyplot as plt
    import numpy as np
    f = np.loadtxt(ffile, usecols=(0,))    

    plt.hist(f, bins=20, color='gray')
    plt.xlabel('$f_{i}$', fontsize=20)
    plt.ylabel('$N_{gc}$', fontsize=20)
    plt.show()
    
    
################################ main program #####################################
gal = raw_input('Type Galaxy number:')

galcatgc = 'N'+gal+'GC.dat'    
galcatpne = 'N'+gal+'PNE.dat'
galcatgc = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatgc
galcatpne = '/home/emilio/MLE/Galaxies/'+gal+'/'+galcatpne

RA, DEC, V, col = read_catalog_gc(galcatgc)
RAp, DECp, Vp = read_catalog_pne(galcatpne) 

path = '/home/emilio/MLE/Galaxies/'+gal+'/'

hdu_list = fits.open(path+'f.fits')
hdu = hdu_list[0]
gal_x = hdu.header['1_XC']
gal_y = hdu.header['1_YC']
gal_x = ast.literal_eval(gal_x)
gal_y = ast.literal_eval(gal_y)


gen_cat(RA, DEC, gal_x, gal_y)          #create 2MASS.GC.cat, the catalog with pixel positions of the tracers.

if os.path.exists(path+"f.fits") == False:
    iraf.imarith(path+"snb.fits", "/", path+"tnb.fits", path+"f.fits")         #create the f.fits, in case it do not exists yet

iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)              #loads some needed iraf libraries
iraf.daophot(_doprint=0)
iraf.ptools(_doprint=0)

#starting photometry

#iraf.phot(path+"snb.fits",path+"2MASS.GC.cat",path+'N'+gal+'B_GC.mag',skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio='yes',sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean='no',rclean=0.,rclip=0.,kclean=0.,mkcenter='no',fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth='no',rgrow=0.,mksky='no',photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert='no',interactive='no',radplots='no',verify='no',update='no',verbose='no', graphics="stdgraph",display="stdimage",icommands="",gcommands="")
#iraf.phot(path+"tnb.fits",path+"2MASS.GC.cat",path+'N'+gal+'T_GC.mag',skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio='yes',sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean='no',rclean=0.,rclip=0.,kclean=0.,mkcenter='no',fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth='no',rgrow=0.,mksky='no',photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert='no',interactive='no',radplots='no',verify='no',update='no',verbose='no', graphics="stdgraph",display="stdimage",icommands="",gcommands="")
iraf.phot(path+"f.fits",path+"2MASS.GC.cat",path+'N'+gal+'F_GC.mag',skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio='yes',sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean='no',rclean=0.,rclip=0.,kclean=0.,mkcenter='no',fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth='no',rgrow=0.,mksky='no',photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert='no',interactive='no',radplots='no',verify='no',update='no',verbose='no', graphics="stdgraph",display="stdimage",icommands="",gcommands="")
 
#dump outputs from iraf in readable format

#iraf.txdump(path+'N'+gal+'B_GC.mag',"XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",'yes', headers='no',paramet='no', Stdout=path+'N'+gal+'B_GC.txt')
#iraf.txdump(path+'N'+gal+'T_GC.mag',"XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",'yes', headers='no',paramet='no', Stdout=path+'N'+gal+'T_GC.txt')
iraf.txdump(path+'N'+gal+'F_GC.mag',"XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",'yes', headers='no',paramet='no', Stdout=path+'N'+gal+'F_GC.txt')

#run ffinder to create the inputs to prob_bulge_bd.f

#ffinder(path+'N'+gal+'B_GC.txt', path+'N'+gal+'T_GC.txt', path+'N'+gal+'F_GC.txt', path+'N'+gal+'_f.GC.dat')
ffinder(path+'N'+gal+'F_GC.txt', path+'N'+gal+'_f.GC.dat')

#plot the fmap if requested:

plot_histf(path+'N'+gal+'_f.GC.dat')

op = raw_input("Plot the fmap? Y/N \n")

if op == "y" or "yes":
    plot_fmap(gal_x, gal_y)
    plot_vmap(gal_x, gal_y)
else: print("Done!")


