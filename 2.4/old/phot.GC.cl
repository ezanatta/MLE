################NGC2768####################################################

## use GALFIT to fit and to create a model disk vs bulge##

##galfit.feedme##

##########################################################################
Note: the coordinate in pixel as obtained from RA and Dec from Vincenzo's catalogue. RA and Dec coordinate in arcsec are the correct one and don't need do be
transfer in pixels because pixel scale 1:1 in 2MASS

############################################################################
Consistency check

displ(image="snb.fits",frame=1,zrange=no,zscale=no,z1=0,z2=1) 
tvmark(frame=1,coords="2MASS.GC.cat",color=205,lab-,num+)


displ(image="tnb.fits",frame=2) 
tvmark(frame=2,coords="2MASS.GC.cat",color=205,lab-,num+)

displ(image="f.fits",frame=3,zrange=no,zscale=no,z1=0,z2=1) 
tvmark(frame=3,coords="2MASS.GC.cat",color=205,lab-,num+)


create 2MASS.GC.cat from prop.f using vprofile.gc.dat

###########################################################################

I checked that they are correctly rotated....

##########################################################################
tnb.fits
snb.fits

############################################################################



## to measure the relative flux use phot (daophot) ##with centerpars=none!!!###

I used sigma=stdv(skys)=0  ===> no bg 
fwhmpsf=2.92 as in galfit psf2.fits (imexam A enclosed)
skyvalue=0 ===> no bg

no centerpars
apertures=1.8 pixels (=1.8 arcsec)

zmag=20


snb.fits & tnb.fits without bg (when I made the model)


!rm N2768B.mag
!rm N2768B.txt
!rm N2768T.mag
!rm N2768T.txt
!rm N2768F.mag
!rm N2768F.txt

!rm N2768B_GC.mag
!rm N2768B_GC.txt
!rm N2768T_GC.mag
!rm N2768T_GC.txt


noao
digiphot
daophot

phot("snb.fits","2MASS.GC.cat","N2768B_GC.mag",skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio=yes,sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean=no,rclean=0.,rclip=0.,kclean=0.,mkcenter=no,fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth=no,rgrow=0.,mksky=no,photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert=no,interactive=no,radplots=no,verify=no,update=no,verbose=no, graphics="stdgraph",display="stdimage",icommands="",gcommands="")

phot("tnb.fits","2MASS.GC.cat","N2768T_GC.mag",skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio=yes,sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean=no,rclean=0.,rclip=0.,kclean=0.,mkcenter=no,fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth=no,rgrow=0.,mksky=no,photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert=no,interactive=no,radplots=no,verify=no,update=no,verbose=no, graphics="stdgraph",display="stdimage",icommands="",gcommands="")

phot("f.fits","2MASS.GC.cat","N2768F_GC.mag",skyfile="",datapars="",scale=1.,fwhmpsf=2.92,emissio=yes,sigma=0.0,datamin=-100, datamax=1e5,noise="poisson",centerpars="",calgorithm="none",cbox=0,cthreshold=0.,minsnratio=1.,cmaxiter=0,maxshift=0,clean=no,rclean=0.,rclip=0.,kclean=0.,mkcenter=no,fitskypars="",salgorithm="constant",annulus=10,dannulus=10,skyvalue=0., smaxiter=10,sloclip=0.,shiclip=0.,snreject=50,sloreject=3.,shireject=3.,khist=3.,binsize=0.1,smooth=no,rgrow=0.,mksky=no,photpars="",weighting="constant",apertures=1.8,zmag=20.,mkapert=no,interactive=no,radplots=no,verify=no,update=no,verbose=no, graphics="stdgraph",display="stdimage",icommands="",gcommands="")
 

txdump("N2768B_GC.mag","XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",yes, headers=no,paramet=no, Stdout="N2768B_GC.txt")
txdump("N2768T_GC.mag","XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",yes, headers=no,paramet=no, Stdout="N2768T_GC.txt")
txdump("N2768F_GC.mag","XCENTER,YCENTER,STDEV,FLUX,SUM,AREA,ID",yes, headers=no,paramet=no, Stdout="N2768F_GC.txt")

###########################################################################

use program ffinder.f to obtain a readable list input for likelihood....


############################################################################
NOTE 193 is not in the image add by hand a value that is reasonabkle (nearby) 0.995.
 
cp N2768_f.GC.dat ../..

#######################################################################
 
