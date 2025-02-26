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

file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored
file_list = sorted(os.listdir(file_folder))  # List the files in the folder

'''
mycmap1 = cm.color_tables.aia_color_table(193 * u.angstrom) # Define colormap to be used as aia 193, set the unit to angstrom
image = f.getdata(f'{file_folder}/{file_list[0]}')  # Get the data from the first image in the file list
imagesize = np.shape(image)  # Get the size of the image

# Plotting the fits iamge using AstroPy
plt.imshow(image, cmap=mycmap1, origin="lower")  # Plot the image
plt.show()  # Show the plot

# Creating an plotting a cropped image using AstroPy
imagecrop = image[2047:4095, 2047:4095]  # Crop the image to the central 2048x2048 pixels
plt.imshow(imagecrop, cmap=mycmap1, origin="lower")  # Plot the cropped image
plt.show()  # Show the plot
'''

# Plotting the fits image as a SunPy Map object
sunpy_map = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic interpolation so the solar north is at the top of the image
fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = sunpy_map_rotated)  # Add a subplot to the figure
sunpy_map_rotated.plot()  # Plot the SunPy map object
plt.title('SunPy Map')  # Set the title of the plot
plt.show()  # Show the plot


# Aligning two SunPy maps and plotting them side by side (alignment durign PLOTTING)
sunpy_map1 = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map1_rotated = sunpy_map1.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic interpolation so the solar north is at the top of the image

sunpy_map2 = sm.Map(f'{file_folder}/{file_list[1]}')  # Create a SunPy map object from the second image
sunpy_map2_rotated = sunpy_map2.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic interpolation so the solar north is at the top of the image


fig = plt.figure(figsize = (12,5))  # Create a new figure

ax1 = fig.add_subplot(121, projection = sunpy_map1_rotated)  # Add a subplot to the figure
sunpy_map1_rotated.plot(axes = ax1)  # Plot the SunPy map object

ax2 = fig.add_subplot(122, projection = sunpy_map2_rotated)  # Add a subplot to the figure
sunpy_map2_rotated.plot(axes = ax2, autoalign = True)  # Plot the SunPy map object aligned with the first image (NOT OVERLAID)

plt.show()  # Show the plot
'''
##########################################################
# Aligning two SunPy maps and storing the aligned objects? (not working)
##########################################################
sunpy_map1 = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map1_rotated = sunpy_map1.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic interpolation so the solar north is at the top of the image

sunpy_map2 = sm.Map(f'{file_folder}/{file_list[1]}')  # Create a SunPy map object from the second image
sunpy_map2_rotated = sunpy_map2.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic interpolation so the solar north is at the top of the image

#sunpy_map1_rotated.plot_settings['norm'] = plt.Normalize(-1500, 1500)
#sunpy_map2_rotated.plot_settings['norm'] = plt.Normalize(-1500, 1500)

# Plot both the images side by side
fig = plt.figure(figsize = (12,5))  # Create a new figure
ax1 = fig.add_subplot(121, projection = sunpy_map1_rotated)  # Add a subplot to the figure
sunpy_map1_rotated.plot(axes = ax1)  # Plot the SunPy map object
ax2 = fig.add_subplot(122, projection = sunpy_map2_rotated)  # Add a subplot to the figure
sunpy_map2_rotated.plot(axes = ax2)
plt.show()  # Show the plot

# Reproject the second image onto the WCS of the first image
out_sunpy_map2_rotated = sunpy_map2_rotated.reproject_interp(sunpy_map1_rotated)  # Reproject the second SunPy map object to the first SunPy map object

# Plot the reprojected image side by side with the first image
fig = plt.figure(figsize = (12,5))  # Create a new figure
ax1 = fig.add_subplot(121, projection = sunpy_map1_rotated)  # Add a subplot to the figure
sunpy_map1_rotated.plot(axes = ax1)  # Plot the first SunPy map object
ax2 = fig.add_subplot(122, projection = out_sunpy_map2_rotated)  # Add a subplot to the figure
out_sunpy_map2_rotated.plot(axes = ax2, titles = 'Reprojected HMI image')  # Plot the reprojected second SunPy map object
plt.show()  # Show the plot

# Plot the images over one another
fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = sunpy_map1_rotated)  # Add a subplot to the figure
sunpy_map1_rotated.plot(axes = ax)  # Plot the first SunPy map object
out_sunpy_map2_rotated.plot(axes = ax, alpha = 0.5)  # Plot the reprojected second SunPy map object over the first SunPy map object with transparency
plt.show()  # Show the plot
'''

# Cropping a map using SkyCoord
top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame = sunpy_map1_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame = sunpy_map1_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
cropped_map = sunpy_map1_rotated.submap(bottom_left, top_right = top_right)  # Crop the SunPy map object

fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = cropped_map)  # Add a subplot to the figure
cropped_map.plot()  # Plot the cropped SunPy map object
plt.show()  # Show the plot


# Using image processing tools to increase image contrast
image_data = cropped_map.data  # Get the data from the cropped SunPy map object
image_data_rescaled = exposure.rescale_intensity(image_data, in_range = (0.1, 0.9))  # Rescale the intensity of the image data

# Create a new SunPy map with the contrast-enhanced data
cropped_map_rescaled = sm.Map(image_data_rescaled, cropped_map.meta)  # Create a new SunPy map object with the rescaled data

fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(projection = cropped_map_rescaled)  # Add a subplot to the figure
cropped_map_rescaled.plot()  # Plot the contrast-enhanced SunPy map object
plt.show()  # Show the plot

# Creating a curvature map from processed image