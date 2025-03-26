import numpy as np
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm
import astropy.units as u
import os
import astropy.io.fits as f


file_folder="NOAA_0-Fe/FITS_files"
file_list=sorted(os.listdir(file_folder))
imagesize=np.shape(f.getdata(f'{file_folder}/{file_list[0]}'))
mycmap1=cm.color_tables.aia_color_table(193*u.angstrom)
image=f.getdata(f'{file_folder}/{file_list[0]}')
plt.imshow(np.rot90(image[:,:],3), cmap=mycmap1, origin="lower")
plt.show()
#plt.imshow(np.rot90(np.log10(image[i,:,:]+5),3), cmap=mycmap1, origin="lower")
#plt.show()


center=[round(imagesize[0]/2), round(imagesize[1]/2)]
halfwidth=2047
plt.imshow(np.rot90(image[center[0]-halfwidth:center[0]+halfwidth,center[1]-halfwidth:center[1]+halfwidth],3), cmap=mycmap1, origin="lower")
plt.show()
    

# Go to helioviewer
# find a time with lots of active regions/spots
# download data
# check the rotation vs helioviewer HMI
