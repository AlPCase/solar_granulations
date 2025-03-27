import os
from src.utils import calc_second_derivative, label_objects, calc_centroids
from src.visualisation import plot_results
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt


def getimagemask(file_index):
    """ Load a FITS file, crop it, and create a binary mask. """

    file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored ############# Should make this a global variable
    file_list = sorted(os.listdir(file_folder))  # List the files in the folder

    # Storing FITS image as SunPy Map object
    sunpy_map = sm.Map(f'{file_folder}/{file_list[file_index]}')  # Create a SunPy map object from the first image
    sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image


    # Cropping a map using SkyCoord
    top_right = SkyCoord(25 * u.arcsec, 25 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
    bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame = sunpy_map_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
    cropped_map = sunpy_map_rotated.submap(bottom_left, top_right = top_right)  # Crop the SunPy map object


    # Create a binary mask? pretty please
    binary_mask = np.zeros(cropped_map.data.shape)  # Create an array of zeros with the same shape as the cropped image

    padded_image = np.pad(cropped_map.data, pad_width=1, mode='constant', constant_values=0)    # Pad the image with NaNs to handle edge cases

    for i in range(1, padded_image.shape[0] - 1):
        for j in range(1, padded_image.shape[0] - 1):
            binary_mask[i-1, j-1] = calc_second_derivative(padded_image, i, j)  # Store the result of calc_second_derivative in the binary mask. Subtract 1 from i and j to account for padding.

    return cropped_map, binary_mask


################### Start of Code ########################

cropped_map, binary_mask = getimagemask(0)

labelled_mask, num_features = label_objects(binary_mask)
centroids = np.array(calc_centroids(labelled_mask))

cropped_map_1, binary_mask_1 = getimagemask(1)

labelled_mask_1, num_features_1 = label_objects(binary_mask_1)
centroids_1 = np.array(calc_centroids(labelled_mask_1))

# Plot the cropped map and binary mask side by side
fig, axs = plt.subplots(2, 3, figsize=(15, 10))  # Create a figure with two rows and three columns of subplots

# Plot the binary maps with centroids of each image, then the centroids placed atop one another
axs[0, 0].imshow(binary_mask, cmap='gray')  # Display the binary mask in the first subplot
axs[0, 0].scatter(centroids[:,1], centroids[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids on the binary mask
axs[0, 0].set_title('OLD Binary Mask with Centroids')  # Set the title of the first subplot

axs[0, 1].imshow(binary_mask_1, cmap='gray')  # Display the binary mask in the second subplot
axs[0, 1].scatter(centroids_1[:,1], centroids_1[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids on the binary mask
axs[0, 1].set_title('NEW Binary Mask with Centroids')  # Set the title of the second subplot

axs[0,2].scatter(centroids[:,1], centroids[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids of the first image
axs[0,2].scatter(centroids_1[:,1], centroids_1[:,0], s=10, c='blue', marker='x', linewidths=1)  # Scatter plot the centroids of the second image

os.makedirs('plots', exist_ok=True)  # Create the 'plots' folder at the root if it doesn't exist
plt.savefig('plots/output_plot_main.png')  # Save the plot in the 'plots' folder with the specified title
