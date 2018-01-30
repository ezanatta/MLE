library(deldir)
library(datautils)

#gal <- paste(readLines('~/MLE/2.3/current_galaxy.dat'), collapse=" ")
gal = 7457
path_to_cat = paste0('~/MLE/Galaxies/',gal, '/')
cat_name = paste0(path_to_cat, '/N',gal,'GC.dat')

cat = read.table(cat_name)

RA = cat$V2
DEC = cat$V3
V = cat$V4
g = cat$V5
i = cat$V6
color = g-i

df = data.frame(RA, DEC)
voronoi = deldir(df$RA, df$DEC)
rbPal <- colorRampPalette(c('blue','red'))
df$Col <- rbPal(10)[as.numeric(cut(V,breaks = 10))]

X11()

plot.deldir(voronoi, fill=df$Col, wlines='tess', wpoints='none', xlab='RA', ylab='DEC', main=paste0('NGC', gal))
dev.copy(jpeg,filename=paste0(path_to_cat,'plots/voronoi', gal,'.jpg'))
