import os
from src.utils import getimagemask, label_objects, trajectories
from src.visualisation import plot_results
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from scipy.spatial import KDTree # For calculating distance between centroids?

# Define max displacement threshold in pixels 
max_displacement = 0.62 # !!!!!!!! Come up with some way for this to be adjusted in km/s and calculated for pixel size (Possible use of FITS header?) !!!!!!!!!!!!

# Number of consective frames to process
num_frames = len(os.listdir("NOAA_0-Fe/FITS_files"))  # Get the number of frames in the folder

# Define list to store centroids for all the fames
all_centroids = []

# Process each frame and calulate the centroids
for frame_index in range(num_frames):
    cropped_map, binary_mask = getimagemask(frame_index)    # Store the FITS image in the file as a cropped sunpy map & create a binary mask
    labelled_mask, num_features = label_objects(binary_mask)    # Label the objects in the binary mask
    centroids = np.array(center_of_mass(labelled_mask, labels=labelled_mask, index=np.unique(labelled_mask)[1:])) # Calculate the centroids of the labelled objects and turn them into a NumPy array
    all_centroids.append(centroids)  # Append the centroids to the list

# PLot the results
fig, axs = plt.subplots(1, 1, figsize=(10, 10))

# Plot trajectories for consecutive frames
for frame_index in range(num_frames - 1):
    centroids_OLD = all_centroids[frame_index]
    centroids_NEW = all_centroids[frame_index + 1]

    # Calculate trajectories between consecutive frames
    x0, y0, x1, y1, matches = trajectories(centroids_OLD, centroids_NEW, max_displacement)  # Calculate the trajectories of the centroids between the two images

    # Plot the centroids
    axs.scatter(centroids_OLD[:, 1], centroids_OLD[:, 0], s=10, c='red', marker='x', label='Frame 1 Centroids')
    axs.scatter(centroids_NEW[:, 1], centroids_NEW[:, 0], s=10, c='blue', marker='x', label='Frame 2 Centroids')

    # Plot the trajectories
    for i in range(len(x0)):
        axs.plot([x0[i], x1[i]], [y0[i], y1[i]], 'k', linewidth=0.5)  # Dashed line for trajectory

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'trajectories_plot.png')
plt.savefig(output_path)
print(f"Plot saved to {output_path}")
