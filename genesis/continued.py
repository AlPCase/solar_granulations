import numpy as np
import matplotlib.pyplot as plt
import os
import sunpy.visualization.colormaps as cm  # Include the SunPy library for use of colormaps
import sunpy.map as sm
import astropy.io.fits as f                 # Include the library for handling FITS files
import astropy.units as u                   # Include the library for units
from astropy.coordinates import SkyCoord    # Include the library for sky coordinates
from reproject import reproject_interp      # Include the library for reprojecting images
from skimage import exposure                # Include the library for image processing
from skimage.filters import sobel
from skimage.feature import canny

file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored
file_list = sorted(os.listdir(file_folder))  # List the files in the folder

# Plotting the fits image as a SunPy Map object
sunpy_map = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image
fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = sunpy_map_rotated)  # Add a subplot to the figure
sunpy_map_rotated.plot()  # Plot the SunPy map object
plt.title('SunPy Map')  # Set the title of the plot
plt.show()  # Show the plot


# Cropping a map using SkyCoord
top_right = SkyCoord(25 * u.arcsec, 25 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
cropped_map = sunpy_map_rotated.submap(bottom_left, top_right = top_right)  # Crop the SunPy map object

fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = cropped_map)  # Add a subplot to the figure
cropped_map.plot()  # Plot the cropped SunPy map object
plt.title('Cropped SunPy Map')  # Set the title of the plot
plt.show()  # Show the plot


# Using image processing tools to increase image contrast
image_data = cropped_map.data  # Get the data from the cropped SunPy map object

# Determine appropriate in_range and out_range for rescaling the image data
p2, p98 = np.percentile(image_data, (2, 98))  # Calculate the 2nd and 98th percentiles of the image data
in_range = (p2, p98)  # Define the in_range for rescaling the image data
out_range = (0, 255)  # Define the out_range for rescaling the image data

image_data_rescaled = exposure.rescale_intensity(image_data, in_range, out_range)  # Rescale the intensity of the image data

# Create a new SunPy map with the contrast-enhanced data
cropped_map_rescaled = sm.Map(image_data_rescaled, cropped_map.meta)  # Create a new SunPy map object with the rescaled data

fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = cropped_map_rescaled)  # Add a subplot to the figure
cropped_map_rescaled.plot()  # Plot the contrast-enhanced SunPy map object
plt.title('Contrast Enhanced SunPy Map')  # Set the title of the plot
plt.show()  # Show the plot


# Create a binary map using curvature criterion - ATTEMPT 1 (UNSUCCESSFUL)
edges = sobel(cropped_map.data)
binary_map = canny(edges)

fig = plt.figure()
plt.imshow(binary_map, cmap = 'gray')
plt.title('Binary Map using Curvature Criterion')
plt.show()