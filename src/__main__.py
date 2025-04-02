import os
from src.utils import calc_second_derivative, label_objects
from src.visualisation import plot_results
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from scipy.spatial import KDTree # For calculating distance between centroids?


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

cropped_map_0, binary_mask_0 = getimagemask(0) # Store the first FITS image in the file as a cropped sunpy map

labelled_mask_0, num_features_0 = label_objects(binary_mask_0)
centroids_0 = np.array(center_of_mass(labelled_mask_0, labels=labelled_mask_0, index=np.unique(labelled_mask_0)[1:])) # Calculate the centroids of the labelled objects and turn them into a NumPy array

cropped_map_1, binary_mask_1 = getimagemask(1)  # Store the second FITS image as a cropped sunpy map

labelled_mask_1, num_features_1 = label_objects(binary_mask_1)
centroids_1 = np.array(center_of_mass(labelled_mask_1, labels=labelled_mask_1, index=np.unique(labelled_mask_1)[1:]))


# Define max displacement threshold in pixels !!!!!!!! This needs to actually be calculated !!!!!!!!!!!!
max_displacement = 2

# Build KDTree from NEW centroids
tree = KDTree(centroids_1[:, ::-1]) # KDTree requires the coordinates to be in (x, y) format, so we reverse the order of the columns

# Query all OLD centroids for nearest neighbours
distances, indices = tree.query(centroids_0[:, ::-1], workers=-1)  # Query the KDTree for the nearest neighbours of each OLD centroid

# Collect all valid matches (i = index in OLD, j = index in NEW)
potential_matches = []
for i, (d,j) in enumerate(zip(distances, indices)):
    if j < len(centroids_1):
        potential_matches.append((i, j, d))

# Sort matches by distance (closest first)
potential_matches.sort(key=lambda x: x[2])

# Greedily select pairs to ensure one-to-one matching
matched_i = []
matched_j = []
matches = []
for match in potential_matches:
    i, j, d = match
    if i not in matched_i and j not in matched_j and d < max_displacement:
        matched_i.append(i)
        matched_j.append(j)
        matches.append((i, j, d))


# Plot the cropped map and binary mask side by side
fig, axs = plt.subplots(2, 3, figsize=(15, 10))  # Create a figure with two rows and three columns of subplots

# Plot the binary maps with centroids of each image, then the centroids placed atop one another
axs[0, 0].imshow(binary_mask_0, cmap='gray')  # Display the binary mask in the first subplot
axs[0, 0].scatter(centroids_0[:,1], centroids_0[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids on the binary mask
axs[0, 0].set_title('OLD Binary Mask with Centroids')  # Set the title of the first subplot

axs[0, 1].imshow(binary_mask_1, cmap='gray')  # Display the binary mask in the second subplot
axs[0, 1].scatter(centroids_1[:,1], centroids_1[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids on the binary mask
axs[0, 1].set_title('NEW Binary Mask with Centroids')  # Set the title of the second subplot

axs[0,2].scatter(centroids_0[:,1], centroids_0[:,0], s=10, c='red', marker='x', linewidths=1)  # Scatter plot the centroids of the first image
axs[0,2].scatter(centroids_1[:,1], centroids_1[:,0], s=10, c='blue', marker='x', linewidths=1)  # Scatter plot the centroids of the second image

# Add lines between matched centroids
for i, j, d in matches:
    x0, y0 = centroids_0[i][1], centroids_0[i][0]  # (y, x) -> (x, y) format
    x1, y1 = centroids_1[j][1], centroids_1[j][0]  # (y, x) -> (x, y) format
    axs[0,2].plot([x0, x1], [y0, y1], 'k', linewidth=1.0)  # Draw a dashed line between the matched centroids


os.makedirs('plots', exist_ok=True)  # Create the 'plots' folder at the root if it doesn't exist
plt.savefig('plots/output_plot_main.png')  # Save the plot in the 'plots' folder with the specified title
