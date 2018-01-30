library(deldir)
library(datautils)

gal <- readline('input galaxy number')
path_to_cat = paste0('~/MLE/Galaxies/',gal, '/')

op <- readline('PNE or GC? ')
op2 <- readline('Bin by Velocity or color? V/c ')

if (op == 'PNE'){
  cat_name = paste0(path_to_cat, '/N',gal,'PNE.dat')
}
if (op == 'GC'){
  cat_name = paste0(path_to_cat, '/N',gal,'GC.dat')
}

cat = read.table(cat_name)

RA = cat$V2
DEC = cat$V3
V = cat$V4
redder = cat$V5
bluer = cat$V6
color = redder-bluer

df = data.frame(RA, DEC)
voronoi = deldir(df$RA, df$DEC)
rbPal <- colorRampPalette(c('blue','red'))

if (op2 == 'V'){
  df$Col <- rbPal(10)[as.numeric(cut(V,breaks = 10))]
}else{
  df$Col <- rbPal(10)[as.numeric(cut(color,breaks = 10))]
}

plot.deldir(voronoi, fill=df$Col, wlines='both', wpoints='none', xlab='RA', ylab='DEC', main=paste0('NGC', gal))
