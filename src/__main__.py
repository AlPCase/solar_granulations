import os
from src.utils import trajectories, velocityField, mraSmoothing, meanDistance
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from scipy.spatial import KDTree
from scipy.interpolate import RBFInterpolator
import pywt


# Define the file folder where the FITS files are stored
file_folder = "NOAA_0-Fe/FITS_files"

# Define max displacement threshold in pixels 
flow_velocity_threshold = 5 # km/s

# Define minimum granule lifetime threshold in seconds (preferably a multiple of cadence)
lifetime_threshold = 135  # seconds

# Define data cadence
cadence = 45  # seconds

# Define the size and location of the region of interest in arcseconds
roi_size = 25  # arcseconds
roi_center = [0, 0] # Bottom-left corner of the region of interest (arcsecond coordinates)


# Number of consective frames to process
#num_frames =  #len(os.listdir(file_folder))  # Get the number of frames in the folder

#for num_frames in range(7, 36):

num_frames = 7

# Total time interval in seconds
temporal_window = num_frames * cadence



trajectories_list = trajectories(num_frames, file_folder, flow_velocity_threshold, roi_size, roi_center)

filtered_trajectories = [
    traj for traj in trajectories_list
    if len(traj["frames"]) * cadence >= lifetime_threshold
]

velocity_field = velocityField(filtered_trajectories, cadence)



##### Wavelet MRA smoothing of the velocity field #####

# Extract positions and velocities from the velocity field
positions = np.array([v["mean_position"] for v in velocity_field])
vx = np.array([v["velocity"][0].value for v in velocity_field])  # Extract x-component of velocity
vy = np.array([v["velocity"][1].value for v in velocity_field])  # Extract y-component of velocity

# Calculate the mean distance between nearest-neighbour points in the velocity field and store to a file for plotting
meanDistance(positions, temporal_window)

# Define grid resolution (must be dyadic?) #! This requires adjustment.
grid_resolution = 32

# Create a dyadic grid for the wavelet transform ????
xi = np.linspace(roi_center[0], roi_center[0] + roi_size*2, grid_resolution)  # X range based on region of interest
yi = np.linspace(roi_center[1], roi_center[1] + roi_size*2, grid_resolution)  # Y range based on region of interest
grid_x, grid_y = np.meshgrid(xi, yi)    # Create a meshgrid for the grid points

# Interpolate the velocities onto the grid using RBFInterpolator
rbf_vx = RBFInterpolator(positions, vx, kernel='thin_plate_spline') # What kernel to use and why?
rbf_vy = RBFInterpolator(positions, vy, kernel='thin_plate_spline')
grid_vx = rbf_vx(np.column_stack([grid_x.ravel(), grid_y.ravel()])).reshape(grid_x.shape)
grid_vy = rbf_vy(np.column_stack([grid_x.ravel(), grid_y.ravel()])).reshape(grid_y.shape)

# Apply wavelet MRA smoothing to the velocity grid
smoothed_vx = mraSmoothing(grid_vx, wavelet='db4', levels=2)
smoothed_vy = mraSmoothing(grid_vy, wavelet='db4', levels=2)



########### Plotting the results ###########
print("\nCalculations complete! Plotting results:")

# Plot the trajectories and velocity field side by side
fig, axs = plt.subplots(2, 2, figsize=(20, 20))  # Create two subplots side by side

# Plot the trajectories
for trajectory in trajectories_list:
    positions = np.array(trajectory["positions"])
    axs[0,0].plot(positions[:, 0], positions[:, 1], marker='x', label=f'Trajectory {trajectory["id"]}')
    axs[0,0].text(positions[-1, 0], positions[-1, 1], str(trajectory["id"]), fontsize=8, color='red')

axs[0,0].set_title("Trajectories")
axs[0,0].set_xlabel("X Position (pixels)")
axs[0,0].set_ylabel("Y Position (pixels)")


# Plot the velocity field
for velocity_data in velocity_field:
    mean_position = velocity_data["mean_position"]
    velocity = velocity_data["velocity"]

    # Plot velocity vector as an arrow
    axs[0,1].quiver(
        mean_position[0], mean_position[1],  # Starting point of the arrow
        velocity[0].value, velocity[1].value,  # Velocity components (x, y)
        angles='xy', scale_units='xy', scale=1, label='Velocity'
    )

axs[0,1].set_title("Velocity Field")
axs[0,1].set_xlabel("X Position (pixels)")
axs[0,1].set_ylabel("Y Position (pixels)")

axs[1,0].quiver(grid_x, grid_y, grid_vx, grid_vy, scale=1, scale_units='xy')
axs[1,0].set_title("Interpolated Velocity Field")

axs[1,1].quiver(grid_x, grid_y, smoothed_vx, smoothed_vy, scale=1, scale_units='xy')
axs[1,1].set_title("Smoothed Velocity Field")


# Add a shared title for the figure
fig.suptitle("Velocity Field", fontsize=16)

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'MRA_velocity_field.png')
plt.savefig(output_path)
print(f"Velocity field plots saved to {output_path}")



# Extract granule lifetimes from trajectories
granule_lifetimes = [len(traj["frames"]) * cadence for traj in trajectories_list]  # Lifetime in seconds

# Define bins based on valid multiples of cadence
max_lifetime = max(granule_lifetimes)  # Maximum lifetime
bins = np.arange(0, max_lifetime + cadence, cadence)  # Bins at intervals of cadence

# Create a histogram of granule lifetimes
plt.figure(figsize=(10, 6))
counts, bin_edges = np.histogram(granule_lifetimes, bins=bins)  # Compute histogram data
plt.step(bin_edges[:-1], counts, where='pre', color='black', linewidth=1.5)  # Step plot

# Add dashed line at y=0
plt.axhline(0, color='black', linestyle='--', linewidth=1)

plt.title("Granule Lifetimes", fontsize=16)
plt.xlabel("Lifetime (seconds)", fontsize=14)
plt.ylabel("Number of Occurrences", fontsize=14)

# Save the plot
os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'granule_lifetimes.png')
plt.savefig(output_path)
print(f"Granule lifetimes plot saved to {output_path}")

print("\n")