import os
from src.utils import getimagemask, label_objects, trajectories
from src.visualisation import plot_results
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from scipy.spatial import KDTree

# Define the file folder where the FITS files are stored
file_folder = "NOAA_0-Fe/FITS_files"

# Define max displacement threshold in pixels 
flow_velocity_threshold = 5 # km/s

# Number of consective frames to process
num_frames = len(os.listdir(file_folder))  # Get the number of frames in the folder

# Initialise trajectories
trajectories_list = []          # List to store trajectories
completed_trajectories = set()  # Set to keep track of completed trajectories
trajectory_id = 0               # Initialise unique ID for each trajectory

# Define list to store centroids for all the fames
all_centroids = []

# Process each frame and calulate the centroids
for frame_index in range(num_frames):
    cropped_map, binary_mask = getimagemask(frame_index, file_folder)        # Store the FITS image in the file as a cropped sunpy map & create a binary mask
    labelled_mask, num_features = label_objects(binary_mask)    # Label the objects in the binary mask
    centroids = np.array(center_of_mass(labelled_mask, labels=labelled_mask, index=np.unique(labelled_mask)[1:])) # Calculate the centroids of the labelled objects and turn them into a NumPy array
    all_centroids.append(centroids)                             # Append the centroids to the list

    if frame_index==0:
        # Initialise the trajectories using the centroids from the first frame
        for trajectory_id, position in enumerate(centroids): #! potential robustness issue Re. trajectory_id?
            trajectories_list.append({
                "id": trajectory_id,        # Unique ID for the trajectory
                "positions": [position],    # List to store the positions of each
                "frames": [frame_index],    # List to store the frames where the trajectory is present
                "index": trajectory_id      # Index of the previous trajectory for linking to next frame
            })
            trajectory_id += 1  # Increment the trajectory ID for the next trajectory
    
    else:
        # Match centroids between OLD and NEW frames

        centroids_OLD = all_centroids[frame_index - 1]  # Centroids from the previous frame
        centroids_NEW = all_centroids[frame_index]      # Centroids from the current frame

        # Calculate the nearest neighbour matched centroids between the two frames
        matches, unmatched_old, unmatched_new = trajectories(centroids_OLD, centroids_NEW, flow_velocity_threshold)  

        # Update the trajectories with the new matches
        matched_trajectories = set()    # Set to keep track of matched trajectories
        for old, new, dist in matches:
    
            trajectory = next(          # Find the trajectory in the list that matches the old centroid
                (
                    t for t in trajectories_list
                    if t["index"] == old
                    and t["id"] not in matched_trajectories     # Prevent adding to the same trajectory multiple times
                    and t["id"] not in completed_trajectories   # Check the trajectory is not "completed"
                    and t["frames"][-1] == frame_index - 1      # Check if the trajectory is from the previous frame (adds robustness against re-adding to "completed" trajectories)
                ),
                None)          # None is returned if no trajectory is found
            
            if trajectory:
                trajectory["positions"].append(centroids_NEW[new])  # Append the new position to the trajectory
                trajectory["frames"].append(frame_index)            # Append the current frame index to the trajectory
                trajectory["index"] = int(new)                      # Update the temp index of the trajectory to that of new centroid.
                matched_trajectories.add(trajectory["id"])          # Mark the trajectory as matched

                print(f"OLD centroid {old}, position {centroids_OLD[old]} matched with NEW centroid {new}, position {centroids_NEW[new]} and added to trajectory ID {trajectory['id']} with distance {dist:.2f}")

        matched_trajectories.clear()  # Clear the set for the next frame
        
        # Create new trajectories from unmatched new centroids
        for Unew in unmatched_new:
            trajectories_list.append({
                "id": trajectory_id,                # Unique ID for the trajectory
                "positions": [centroids_NEW[Unew]], # List to store the positions of each centroid
                "frames": [frame_index],            # List to store the frames where the trajectory is present
                "index": Unew                       # Index of the previous trajectory for linking to next frame
            })
            trajectory_id += 1              # Increment the trajectory ID for the next trajectory to be created
        
        # "end" the trajectories corresponding unmatched old centroids
        for Uold in unmatched_old:
            trajectory = next(
                (
                    t for t in trajectories_list
                    if t["index"] == Uold           # Find the trajectory in the list that matches the old centroid
                ),
                None)        # None is returned if no trajectory is found
            
            if trajectory:
                trajectory["index"] = -1                      # Mark the trajectory as ended by setting the index to -1 (no new centroid positions will be added)
                completed_trajectories.add(trajectory["id"])  # Mark the trajectory as completed

        
        print("\n")

print(f"Trajectories: {trajectories_list}")  # Print the list of trajectories


# Plot the trajectories
fig, axs = plt.subplots(1, 1, figsize=(10, 10))

for trajectory in trajectories_list:
    positions = np.array(trajectory["positions"])
    axs.plot(positions[:, 0], positions[:, 1], marker='x', label=f'Trajectory {trajectory["id"]}')
    axs.text(positions[-1, 0], positions[-1, 1], str(trajectory["id"]), fontsize=8, color='red')

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'trajectories_plot.png')
plt.savefig(output_path)
print(f"Plot saved to {output_path}")
