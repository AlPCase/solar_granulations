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
flow_velocity_threshold = 5 # km/s

# Number of consective frames to process
num_frames = len(os.listdir("NOAA_0-Fe/FITS_files"))  # Get the number of frames in the folder

# Initialise trajectories
trajectories_list = [] # List to store trajectories
all_trajectories = {} # Dictionary to store active trajectories
trajectory_id = 0 # Unique ID for each trajectory

# Define list to store centroids for all the fames
all_centroids = []

# Process each frame and calulate the centroids
for frame_index in range(num_frames):
    cropped_map, binary_mask = getimagemask(frame_index)    # Store the FITS image in the file as a cropped sunpy map & create a binary mask
    labelled_mask, num_features = label_objects(binary_mask)    # Label the objects in the binary mask
    centroids = np.array(center_of_mass(labelled_mask, labels=labelled_mask, index=np.unique(labelled_mask)[1:])) # Calculate the centroids of the labelled objects and turn them into a NumPy array
    all_centroids.append(centroids)  # Append the centroids to the list

    if frame_index==0:
        # Initialise the trajectories using the centroids from the first frame
        for trajectory_id, position in enumerate(centroids):
            trajectories_list.append({
                "id": trajectory_id,        # Unique ID for the trajectory
                "positions": [position],    # List to store the positions of each
                "frames": [frame_index],    # List to store the frames where the trajectory is present
                "index": trajectory_id      # Index of the previous trajectory for linking to next frame
            })
            all_trajectories[trajectory_id] = trajectories_list[-1] # Store the trajectory in the dictionary
            trajectory_id += 1  # Increment the trajectory ID for the next trajectory
    
    else:
        # Match centroid between OLD and NEW frames
        centroids_OLD = all_centroids[frame_index - 1]
        centroids_NEW = all_centroids[frame_index]
        x0, y0, x1, y1, matches = trajectories(centroids_OLD, centroids_NEW, flow_velocity_threshold)  # Calculate the nearest neighbour matched centroids between the two frames

        # Update the trajectories with the new matches
        for i, (index_OLD, index_NEW, _) in enumerate(matches):
            # Find the trajectory in all_trajectories where "index" matches index_OLD
            trajectory = next((traj for traj in all_trajectories.values() if traj["index"] == index_OLD), None)
            
            trajectory["positions"].append((x1[index_OLD], y1[index_OLD]))  # Add the new centroid position
            trajectory["frames"].append(frame_index)  # Add the current frame index
            trajectory["index"] = index_NEW  # Update the index to the new centroid
            


        # If the granule in the PREVIOUS frame doesn't match any granule in the CURRENT frame,
        # count the trajectory as ENDED? (probably not, I think it will be fine if simply no new matches are added.)

        # If the granule in the CURRENT frame doesn't match any granule in the PREVIOUS frame,
        # start a new trajectory with a unique ID.

print(all_trajectories)

# PLot the results
fig, axs = plt.subplots(1, 1, figsize=(10, 10))

# Plot trajectories for consecutive frames
for frame_index in range(num_frames - 1):
    centroids_OLD = all_centroids[frame_index]
    centroids_NEW = all_centroids[frame_index + 1]

    # Calculate trajectories between consecutive frames
    x0, y0, x1, y1, matches = trajectories(centroids_OLD, centroids_NEW, flow_velocity_threshold)  # Calculate the trajectories of the centroids between the two images

    # Plot the centroids
    axs.scatter(centroids_OLD[:, 1], centroids_OLD[:, 0], s=10, c='red', marker='x', label='Frame 1 Centroids')
    axs.scatter(centroids_NEW[:, 1], centroids_NEW[:, 0], s=10, c='blue', marker='x', label='Frame 2 Centroids')

    # Plot the trajectories
    for i in range(len(x0)):
        axs.plot([x0[i], x1[i]], [y0[i], y1[i]], 'k', linewidth=0.5)  # Plot the trajectory lines between the centroids

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'trajectories_plot.png')
plt.savefig(output_path)
print(f"Plot saved to {output_path}")
