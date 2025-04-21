import os
import sunpy.map as sm
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import label, center_of_mass
from scipy.spatial import KDTree  # For determining nearest neighbour centroids
import matplotlib.pyplot as plt

def trajectories(centroids_OLD, centroids_NEW):
    """ Match centroids between consecutive frames. """


    # Create a KDTree from NEW centroids
    tree = KDTree(centroids_NEW[:, :2]) # KDTree requires 

    # Query all OLD centroids for nearest neighbours in NEW centroids
    distances, indices = tree.query(centroids_OLD[:, :2])  # Find the nearest neighbour for each OLD centroid

    # Collect all potential matches with their distances
    potential_matches = []
    for old, (dist, new) in enumerate(zip(distances, indices)):
        potential_matches.append((old, new, dist)) # old refers to the index in OLD, new refers to the index in NEW



    # Sort matches by distance (closest first)
    potential_matches.sort(key=lambda x: x[2])

    matched_old = []  # List to keep track of matched OLD indices
    matched_new = []  # List to keep track of matched NEW indices
    matches = []  # List to store the final matches
    for match in potential_matches:
        old, new, dist = match
        if old not in matched_old and new not in matched_new:
            matched_old.append(old)
            matched_new.append(new)
            matches.append((old, new, dist))

    # Create a list of unmatched OLD and NEW centroids
    unmatched_old = [i for i in range(len(centroids_OLD)) if i not in matched_old]
    unmatched_new = [i for i in range(len(centroids_NEW)) if i not in matched_new]
    
    return matches, unmatched_old, unmatched_new




######### START OF CODE ###########

# Define generic dataset of centroids for frames 1 and 2
centroids_0 = np.array([[1, 2], [3, 4], [5, 6]])
centroids_1 = np.array([[3.1, 4.1], [1.1, 2.1], [5.1, 6.1], [7.1, 8.1]])
centroids_2 = np.array([[1.2, 2.2], [3.2, 4.2], [5.2, 6.2]])

# Define list to store centroids for all the frames
ID_centroids = []

# Append the ID_centroids list with the centroids from each frame, assigning unique IDs to each centroid
centroids_0_with_ID = [(i, centroid) for i, centroid in enumerate(centroids_0)]
ID_centroids.append(centroids_0_with_ID)  # Append the centroids from the first frame to the list

centroids_1_with_ID = [(i, centroid) for i, centroid in enumerate(centroids_1)]
ID_centroids.append(centroids_1_with_ID)  # Append the centroids from the second frame to the list

centroids_2_with_ID = [(i, centroid) for i, centroid in enumerate(centroids_2)]
ID_centroids.append(centroids_2_with_ID)  # Append the centroids from the third frame to the list

# Define lists for storing trajectories
trajectories_list = []  # List to store unique trajectory data
trajectory_ID = 0       # Unique ID for each trajectory, starting at 0.


# Print the centroids with IDs for each frame.
for frame_index, ID_and_centroid in enumerate(ID_centroids): # This enumeration is only relevant in this test code I think.
    print(f"Frame {frame_index}:")

    
    # Initialise trajectories using the centroids from the first frame
    if frame_index == 0:
        for i in ID_and_centroid:
            print(f"Centroid ID: {i[0]}, Centroid Coordinates: {i[1]}")
            trajectories_list.append({
                "id": trajectory_ID,        # Unique ID for the trajectory
                "positions": [i[1]],        # List to store the positions of each centroid
                "frames": [frame_index],    # List to store the frames where the trajectory is present
                "index": i[0]      # Index of the previous trajectory for linking to next frame
            })
            trajectory_ID += 1  # Increment the trajectory ID for the next trajectory

    else:
        for i in ID_and_centroid:
            print(f"Centroid ID: {i[0]}, Centroid Coordinates: {i[1]}")

        # Extract the centroids from the previous and current frames
        centroids_old = np.array([centroid for _, centroid in ID_centroids[frame_index - 1]])
        centroids_new = np.array([centroid for _, centroid in ID_centroids[frame_index]])

        matches, unmatched_old, unmatched_new = trajectories(centroids_old, centroids_new)
        # For the centroid in the NEW frame with matching ID to the OLD frame, update the trajectory
        
        matched_trajectories = set() # Set to keep track of matched trajectories

        # Update trajectories based on matches
        for old, new, dist in matches:

            print(f"OLD centroid {old}, position {centroids_old[old]} matched with NEW centroid {new}, position {centroids_new[new]} with distance {dist:.2f}")

            # Find the trajectory corresponding to the OLD centroid
            trajectory = next(
                        (
                            traj for traj in trajectories_list 
                            if traj["index"] == old and traj["id"] not in matched_trajectories
                        ),
                        None)
            
            if trajectory:
                trajectory["positions"].append(centroids_new[new])  # Add the new centroid position
                trajectory["frames"].append(frame_index)            # Add the current frame index
                trajectory["index"] = int(new)                      # Update the index to the NEW centroid index
                matched_trajectories.add(trajectory["id"])          # Mark the trajectory as matched

                print(f"Updated trajectory ID {trajectory['id']} with new position {centroids_new[new]} and frame {frame_index}")

        matched_trajectories.clear()  # Clear the set for the next frame
        

        # Create new trajectories from unmatched new centroids
        for Unew in unmatched_new:
            trajectories_list.append({
                "id": trajectory_ID,        # Unique ID for the trajectory
                "positions": [centroids_new[Unew]],  # List to store the positions of each centroid
                "frames": [frame_index],    # List to store the frames where the trajectory is present
                "index": Unew                # Index of the previous trajectory for linking to next frame
            })
            trajectory_ID += 1

            print(f"Created new trajectory ID {trajectory_ID - 1} with position {centroids_new[Unew]} and frame {frame_index}")


        for Uold in unmatched_old:
            # Check if the OLD centroid is part of any trajectory
            trajectory = next((traj for traj in trajectories_list if traj["index"] == Uold), None)
            if trajectory:
                trajectory["index"] = None # Set the temp index to None to indicate that the trajectory is completed and prevent future updates.

                print(f"Completed trajectory ID {trajectory['id']} with position {centroids_old[Uold]} and frame {frame_index - 1}")

    print("\n")

print("Stored Trajectories:")
for trajectory in trajectories_list:
    print(trajectory)


# Plot the trajectories for the first two frames
fig, axs = plt.subplots(1, 1, figsize=(10, 10))

for trajectory in trajectories_list:
    positions = np.array(trajectory["positions"])
    axs.plot(positions[:, 0], positions[:, 1], marker='x', label=f'Trajectory {trajectory["id"]}')
    axs.text(positions[0, 0], positions[0, 1], f'Start {trajectory["id"]}', fontsize=8, color='red')


os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'test.png')
plt.savefig(output_path)
print(f"Plot saved to {output_path}")

