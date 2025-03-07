import numpy as np
import matplotlib.pyplot as plt
import os
import sunpy.visualization.colormaps as cm  # Include the SunPy library for use of colormaps
import sunpy.map as sm
import astropy.io.fits as f                 # Include the library for handling FITS files
import astropy.units as u                   # Include the library for units
from astropy.coordinates import SkyCoord    # Include the library for sky coordinates
from reproject import reproject_interp      # Include the library for reprojecting images


def calc_second_derivative(image, i, j):
    
    # Define a 3x3 kernel around the pixel of interest that includes adjacent pixels
    kernel_image = ([image[i-1, j-1], image[i-1, j], image[i-1, j+1]],
                    [image[i, j-1], image[i, j], image[i, j+1]],
                    [image[i+1, j-1], image[i+1, j], image[i+1, j+1]])

    # Define kernels that approximate second derivative using laplacian operator
    kernel_ns = np.array([[0, 1, 0], # North -> South
                        [0, -2, 0],
                        [0, 1, 0]])

    kernel_ew = np.array([[0, 0, 0], # East -> West
                        [1, -2, 1],
                        [0, 0, 0]])

    kernel_nesw = np.array([[0, 0, 1],      # North East -> South West
                            [0, -2, 0], 
                            [1, 0, 0]])

    kernel_nwse = np.array([[1, 0, 0],      # North West -> South East
                            [0, -2, 0], 
                            [0, 0, 1]])

    # Calculate second derivatives
    sd_ns = np.sum(kernel_image * kernel_ns)
    sd_ew = np.sum(kernel_image * kernel_ew)
    sd_nesw = np.sum(kernel_image * kernel_nesw)
    sd_nwse = np.sum(kernel_image * kernel_nwse)

    # If the signs of all the second derivatives are the same, set this pixel as an "object pixel"
    if np.sign(sd_ns) == np.sign(sd_ew) == np.sign(sd_nesw) == np.sign(sd_nwse):
        if np.sign(sd_ns) == -1:    # If the sign of the second derivative is negative, the object is called "bright".
            return 1
        else:
            return 0    # otherwise, the object is called "dark".


file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored
file_list = sorted(os.listdir(file_folder))  # List the files in the folder


# Storing FITS image as SunPy Map object
sunpy_map = sm.Map(f'{file_folder}/{file_list[0]}')  # Create a SunPy map object from the first image
sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image


# Cropping a map using SkyCoord
top_right = SkyCoord(25 * u.arcsec, 25 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
cropped_map = sunpy_map_rotated.submap(bottom_left, top_right = top_right)  # Crop the SunPy map object


# Create a binary mask? pretty please
binary_mask = np.zeros(cropped_map.data.shape)  # Create an array of zeros with the same shape as the cropped image

padded_image = np.pad(cropped_map.data, pad_width=1, mode='constant', constant_values=0)    # Pad the image with zeros to handle edge cases

for i in range(1, padded_image.shape[0] - 1):
    for j in range(1, padded_image.shape[0] - 1):
        binary_mask[i-1, j-1] = calc_second_derivative(padded_image, i, j)  # Store the result of calc_second_derivative in the binary mask. Subtract 1 from i and j to account for padding.


# Plot the cropped map and binary mask side by side
fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # Create a figure with two subplots
axs[0].imshow(cropped_map.data, cmap = 'gray')  # Display the cropped map in the first subplot
axs[0].set_title('Cropped Map')  # Set the title of the first subplot
axs[1].imshow(binary_mask, cmap = 'gray')  # Display the binary mask in the second subplot
axs[1].set_title('Binary Mask')  # Set the title of the second subplot
axs[2].imshow(binary_mask, cmap ='gray')
axs[2].imshow(cropped_map.data, cmap = 'gray', alpha=0.5)  # Display the cropped map in the third subplot with 50% opacity
axs[2].set_title('Binary Mask Overlay')  # Set the title of the third subplot
plt.show()