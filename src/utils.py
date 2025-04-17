import os
import sunpy.map as sm
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import label, center_of_mass
from scipy.spatial import KDTree  # For determining nearest neighbour centroids


def calc_second_derivative(image, i, j):
    """Calculate the second derivative for a given pixel."""

    # Define a 3x3 kernel around the pixel of interest that includes adjacent pixels
    kernel_image = ([image[i-1, j-1], image[i-1, j], image[i-1, j+1]],
                    [image[i, j-1], image[i, j], image[i, j+1]],
                    [image[i+1, j-1], image[i+1, j], image[i+1, j+1]])
    
    # Define kernels that approximate second derivative using laplacian operator
    kernel_ns = np.array([[0, 1, 0],
                          [0, -2, 0], 
                          [0, 1, 0]])  # North-South
    
    kernel_ew = np.array([[0, 0, 0],
                          [1, -2, 1],
                          [0, 0, 0]])  # East-West
    
    kernel_nesw = np.array([[0, 0, 1], 
                            [0, -2, 0],
                            [1, 0, 0]])  # NE-SW
    
    kernel_nwse = np.array([[1, 0, 0],
                            [0, -2, 0],
                            [0, 0, 1]])  # NW-SE
    
    # Calculate second derivatives
    sd_ns = np.sum(kernel_image * kernel_ns)
    sd_ew = np.sum(kernel_image * kernel_ew)
    sd_nesw = np.sum(kernel_image * kernel_nesw)
    sd_nwse = np.sum(kernel_image * kernel_nwse)

    # If the signs of all the second derivatives are the same, set this pixel as an "object pixel"
    if np.sign(sd_ns) == np.sign(sd_ew) == np.sign(sd_nesw) == np.sign(sd_nwse) == -1:
        return 1  # Bright object pixel (belongs to a granule)
    else:
        return 0  # Dark object pixel


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



def label_objects(binary_mask):
    """Label connected components in the binary mask."""

    labelled_mask = np.zeros(binary_mask.shape)  # Create an array of zeros with the same shape as the binary mask

    s = np.ones((3, 3), dtype=np.int32)  # Structuring element, 3x3 square
    labelled_mask, num_features = label(binary_mask, structure=s) # Create an array of zeros with the same shape as the binary mask

    return labelled_mask, num_features


def trajectories(centroids_OLD, centroids_NEW, max_displacement):
    """Calculate the trajectories of centroids between two images."""
    
    # Build KDTree from NEW centroids
    tree = KDTree(centroids_NEW[:, ::-1]) # KDTree requires the coordinates to be in (x, y) format, so we reverse the order of the columns

    # Query all OLD centroids for nearest neighbours
    distances, indices = tree.query(centroids_OLD[:, ::-1])  # Query the KDTree for the nearest neighbours of each OLD centroid

    # Collect all valid matches (i = index in OLD, j = index in NEW)
    potential_matches = []
    for i, (d, j) in enumerate(zip(distances, indices)):
        if d <= max_displacement:  # Only include matches within the threshold
            potential_matches.append((i, j, d))

    # Sort matches by distance (x[2]) (closest first)
    potential_matches.sort(key=lambda x: x[2])

    # Ensure one-to-one matching by checking that the centroids are not already matched, choosing the matches with the shortest distance first.
    matched_i = []
    matched_j = []
    matches = []
    for match in potential_matches:
        i, j, d = match # Potentially redundant line? Could just write for i, j, d in potential_matches: ????????????
        if i not in matched_i and j not in matched_j and d < max_displacement:
            matched_i.append(i)
            matched_j.append(j)
            matches.append((i, j, d))

    # Extract coordinates of matched centroids
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    for i, j, d in matches:
        y0.append(centroids_OLD[i][0])  # y-coordinate of centroids_0
        x0.append(centroids_OLD[i][1])  # x-coordinate of centroids_0
        y1.append(centroids_NEW[j][0])  # y-coordinate of centroids_1
        x1.append(centroids_NEW[j][1])  # x-coordinate of centroids_1

    # Convert to NumPy arrays
    x0 = np.array(x0)
    y0 = np.array(y0)
    x1 = np.array(x1)
    y1 = np.array(y1)


    return x0, y0, x1, y1, matches