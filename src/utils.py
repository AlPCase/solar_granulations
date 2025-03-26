import numpy as np
from scipy.ndimage import label, center_of_mass

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

def label_objects(binary_mask):
    """Label connected components in the binary mask."""

    s = np.ones((3, 3), dtype=np.int32)  # Structuring element, 3x3 square
    labelled_mask, num_features = np.zeros(binary_mask.shape) # Create an array of zeros with the same shape as the binary mask

    return labelled_mask, num_features

def calc_centroids(labelled_mask):
    """Calculate centroids of labeled objects."""
    centroids = center_of_mass(labelled_mask, labels=labelled_mask, index=np.unique(labelled_mask)[1:])
    
    return np.array(centroids)