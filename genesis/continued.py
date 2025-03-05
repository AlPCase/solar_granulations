import numpy as np
import matplotlib.pyplot as plt
import os
import sunpy.visualization.colormaps as cm  # Include the SunPy library for use of colormaps
import sunpy.map as sm
import astropy.io.fits as f                 # Include the library for handling FITS files
import astropy.units as u                   # Include the library for units
from astropy.coordinates import SkyCoord    # Include the library for sky coordinates
from reproject import reproject_interp      # Include the library for reprojecting images


file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored
file_list = sorted(os.listdir(file_folder))  # List the files in the folder


# Storing FITS image as SunPy Map object
sunpy_map = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image

# Cropping a map using SkyCoord
top_right = SkyCoord(25 * u.arcsec, 25 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
cropped_map = sunpy_map_rotated.submap(bottom_left, top_right = top_right)  # Crop the SunPy map object


fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = cropped_map)  # Add a subplot to the figure
cropped_map.plot()  # Plot the cropped SunPy map object
plt.title('Cropped SunPy Map')  # Set the title of the plot
plt.show()  # Show the plot
# Add line to leave the plot open and continue code


# Create a new image by stepping through each pixel of the cropped map and creating a new np array
new_image = np.zeros(cropped_map.data.shape)  # Create a new np array of zeros with the same shape as the cropped map

for i in range(cropped_map.data.shape[0]):  # Loop over the rows of the cropped map
    for j in range(cropped_map.data.shape[1]):  # Loop over the columns of the cropped map
        new_image[i, j] = cropped_map.data[i, j]  # Set the value of the new image at the current pixel to the value of the cropped map at the same pixel

fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(111)
ax.imshow(new_image, cmap = 'gray')
plt.title('New Image from Array Indexing')
plt.show()  # Show the plot

