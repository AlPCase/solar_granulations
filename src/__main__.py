import os
from src.utils import calc_second_derivative, label_objects, calc_centroids
from src.visualisation import plot_results
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

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
centroids = calc_centroids(labelled_mask)

plot_results(cropped_map, binary_mask, centroids)
