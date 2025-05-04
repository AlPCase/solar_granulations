import os
from src.utils import getimagemask, label_objects, matching, trajectories
import sunpy.map as sm
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from scipy.spatial import KDTree
from scipy.interpolate import RBFInterpolator
import pywt

def wavelet_mra_smooth(grid, wavelet='db4', levels=4):
    """ Apply wavelet MRA smoothing to a 2D grid of velocities. """
    # Decompose the grid into wavelet coefficients
    coeffs = pywt.wavedec2(grid, wavelet, level=levels)

    # Keep only the approximation coefficients at the coarsest scale
    coeffs_smooth = [coeffs[0]]
    for i in range(1, len(coeffs)):
        coeffs_smooth.append((
            (np.zeros_like(coeffs[i][0])),
            (np.zeros_like(coeffs[i][1])),
            (np.zeros_like(coeffs[i][2]))
        ))

    return pywt.waverec2(coeffs_smooth, wavelet)


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
num_frames = 17 #len(os.listdir(file_folder))  # Get the number of frames in the folder

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
            V = delta_X / delta_t #displacement * u.pix #/ delta_t.value * u.pix/u.s #/ delta_t #! madness
        else:
            V = np.array([0, 0]) * u.pix #u.km / u.s  # Assign zero velocity if no time difference

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


## Wavelet MRA smoothing of the velocity field ##

# Extract positions and velocities from the velocity field
positions = np.array([v["mean_position"] for v in velocity_field])
vx = np.array([v["velocity"][0].value for v in velocity_field])  # Extract x-component of velocity
vy = np.array([v["velocity"][1].value for v in velocity_field])  # Extract y-component of velocity

#print(f"Positions: {positions}")

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
smoothed_vx = wavelet_mra_smooth(grid_vx, wavelet='db4', levels=2)
smoothed_vy = wavelet_mra_smooth(grid_vy, wavelet='db4', levels=2)



########### Plotting the results ###########

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
print(f"Plot saved to {output_path}")



# Extract granule lifetimes from trajectories
granule_lifetimes = [len(traj["frames"]) * dt for traj in trajectories_list]  # Lifetime in seconds

# Define bins based on valid multiples of dt
max_lifetime = max(granule_lifetimes)  # Maximum lifetime
bins = np.arange(0, max_lifetime + dt, dt)  # Bins at intervals of dt

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
