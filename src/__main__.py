import os
from src.utils import getimagemask, label_objects, matching, trajectories
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

# Data time interval
dt = 45  # seconds

# Define max displacement threshold in pixels 
flow_velocity_threshold = 5 # km/s

# Define the size and location of the region of interest in arcseconds
roi_size = 25  # arcseconds
roi_center = [0, 0] # Bottom-left corner of the region of interest (arcsecond coordinates)


# Number of consective frames to process
num_frames = len(os.listdir(file_folder))  # Get the number of frames in the folder

trajectories_list = trajectories(num_frames, file_folder, flow_velocity_threshold, roi_size, roi_center)

velocity_field = []
for traj in trajectories_list:
    positions = np.array(traj["positions"])
    frames = np.array(traj["frames"])


    if len(frames) > 1:
        # Calculate the mean velocity associated with each trajectory
        # V_k = (X_n2 - X_n1) / (t_n2 - t_n1)
        displacement = positions[-1] - positions[0] # Calculate displacement in pixels as the difference between the final and first positions
        delta_X = displacement * (725/2) * u.km  # Convert to km

        # Calculate number of frames and convert to seconds
        active_frames = len(frames) - 1
        delta_t = active_frames * dt * u.s  # Converting to seconds

        if delta_t > 0:
            V = delta_X / delta_t
        else:
            V = np.array([0, 0]) * u.km / u.s  # Assign zero velocity if no time difference

        # Calculate the mean position of each trajectory
        # X_k = 1/(n1 - n2 + 1) * sum(X_n1, X_n2, ..., X_nk)
        mean_position = np.mean(positions, axis=0)  # Mean position of the trajectory

        # Combine the mean velocity and mean position of each trajectory
        # into set {V_k, X_k} for each trajectory k
        velocity_field.append({
            "velocity": V,
            "mean_position": mean_position,
            "trajectory_id": traj["id"]
        })
    
    else:
        continue

# Plot the trajectories and velocity field side by side
fig, axs = plt.subplots(1, 2, figsize=(20, 10))  # Create two subplots side by side

# Plot the trajectories
for trajectory in trajectories_list:
    positions = np.array(trajectory["positions"])
    axs[0].plot(positions[:, 0], positions[:, 1], marker='x', label=f'Trajectory {trajectory["id"]}')
    axs[0].text(positions[-1, 0], positions[-1, 1], str(trajectory["id"]), fontsize=8, color='red')

axs[0].set_title("Trajectories")
axs[0].set_xlabel("X Position (pixels)")
axs[0].set_ylabel("Y Position (pixels)")


# Plot the velocity field
for velocity_data in velocity_field:
    mean_position = velocity_data["mean_position"]
    velocity = velocity_data["velocity"]

    # Plot velocity vector as an arrow
    axs[1].quiver(
        mean_position[0], mean_position[1],  # Starting point of the arrow
        velocity[0].value, velocity[1].value,  # Velocity components (x, y)
        angles='xy', scale_units='xy', scale=1, label='Velocity'
    )

axs[1].set_title("Velocity Field")
axs[1].set_xlabel("X Position (pixels)")
axs[1].set_ylabel("Y Position (pixels)")

# Add a shared title for the figure
fig.suptitle("Unresolved Velocity Field", fontsize=16)

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'unresolved_velocity_field.png')
plt.savefig(output_path)
print(f"Plot saved to {output_path}")

