import wget
import os

def down_index():
    path = 'http://data.astrometry.net/4200/index-42'
    fits = '.fits'    
       
    for i in range(0,7): 
        for j in range(0, 48):
            if j < 10:
                filename = '0'+str(i)+'-0'+str(j)
            else:
                filename = '0'+str(i)+'-'+str(j)
            url = path+filename+fits
            if os.path.exists('/usr/local/astrometry/data/index-42'+filename+fits) == False:
                print url
                index = wget.download(url, out='/usr/local/astrometry/data')
            
down_index()