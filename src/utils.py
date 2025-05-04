import os
import sunpy.map as sm
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import label, center_of_mass
from scipy.spatial import KDTree


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


def getimagemask(file_index, file_folder, roi_size, roi_center):
    """ Load a FITS file, crop it, and create a binary mask. """

    file_list = sorted(os.listdir(file_folder))  # List the files in the folder

    # Storing FITS image as SunPy Map object
    sunpy_map = sm.Map(f'{file_folder}/{file_list[file_index]}')  # Create a SunPy map object from the first image
    sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image


    # Cropping a map using SkyCoord
    top_right = SkyCoord(roi_size*u.arcsec, roi_size*u.arcsec, frame=sunpy_map_rotated.coordinate_frame)  # Define the top right corner of the cropped image as a SkyCoord object
    bottom_left = SkyCoord(roi_center[0]*u.arcsec, roi_center[1]*u.arcsec, frame=sunpy_map_rotated.coordinate_frame)  # Define the bottom left corner of the cropped image as a SkyCoord object
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


def matching(centroids_OLD, centroids_NEW, flow_velocity_threshold):
    """Calculate the trajectories of centroids between two images."""

    # Convert km/s to pixels/frame
    max_displacement = (flow_velocity_threshold * 45) / (725 * 0.5)  # 45s temporal relution, 1" (effective resolution) is 725km, each pixel is 0.5"
    
    # Create a KDTree from NEW centroids
    tree = KDTree(centroids_NEW[:, :2]) # KDTree requires 

    # Query all OLD centroids for nearest neighbours in NEW centroids
    distances, indices = tree.query(centroids_OLD[:, :2])  # Find the nearest neighbour for each OLD centroid

    # Collect all potential matches with their distances
    potential_matches = []
    for old, (dist, new) in enumerate(zip(distances, indices)):
        if dist < max_displacement:
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
            matches.append((old, int(new), dist))

    # Create a list of unmatched OLD and NEW centroids
    unmatched_old = [i for i in range(len(centroids_OLD)) if i not in matched_old]
    unmatched_new = [i for i in range(len(centroids_NEW)) if i not in matched_new]
    
    return matches, unmatched_old, unmatched_new


def trajectories(num_frames, file_folder, flow_velocity_threshold, roi_size, roi_center):
    """ Calculate the trajectories of centroids in a series of images. """

    # Initialise trajectories
    trajectories_list = []          # List to store trajectories
    completed_trajectories = set()  # Set to keep track of completed trajectories
    trajectory_id = 0               # Initialise unique ID for each trajectory

    # Define list to store centroids for all the fames
    all_centroids = []

    # Process each frame and calulate the centroids
    for frame_index in range(num_frames):
        cropped_map, binary_mask = getimagemask(frame_index, file_folder, roi_size, roi_center)        # Store the FITS image in the file as a cropped sunpy map & create a binary mask
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
            matches, unmatched_old, unmatched_new = matching(centroids_OLD, centroids_NEW, flow_velocity_threshold)  

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

                    #print(f"OLD centroid {old}, position {centroids_OLD[old]} matched with NEW centroid {new}, position {centroids_NEW[new]} and added to trajectory ID {trajectory['id']} with distance {dist:.2f}")

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

    return trajectories_list