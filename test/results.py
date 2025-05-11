import os
from utils import trajectories, velocityField
import numpy as np
import pandas as pd


# Define the file folder where the FITS files are stored
file_folder = r"C:/Users/AlexC/source/repos/solar_granulations/NOAA_0-Fe/FITS_files"

# Define Parameters
flow_velocity_threshold = 5 # km/s
lifetime_threshold = 135  # seconds, granule lifetime threshold
cadence = 45  # seconds
num_frames = 36 #len(os.listdir(file_folder))  # Get the number of frames in the folder
roi_size = 50  # arcseconds
roi_centers = [[0, y] for y in range(-1000, 1001, 50)]

roi_center = [-400, -25]

# Total time interval in seconds
temporal_window = num_frames * cadence

#for roi_center in roi_centers:

#roiIndex = roi_centers.index(roi_center) + 1
#print(f"Processing ROI {roiIndex} of {len(roi_centers)} Location: {roi_center}")

# Determine the trajectories across the dataset
trajectories_list = trajectories(num_frames, file_folder, flow_velocity_threshold, roi_size, roi_center)

"""# Filter the trajectories based on the lifetime threshold
filtered_trajectories = [
    traj for traj in trajectories_list
    if len(traj["frames"]) * cadence >= lifetime_threshold
]"""

velocity_field = velocityField(trajectories_list, cadence)

# Extract positions and velocities from the velocity field
positions = np.array([v["mean_position"] for v in velocity_field])
vx = np.array([v["velocity"][0].value for v in velocity_field])  # Extract x-component of velocity
vy = np.array([v["velocity"][1].value for v in velocity_field])  # Extract y-component of velocity

# Calculate the mean velocity in the x and y directions 
meanvx = np.mean(vx)
meanvy = np.mean(vy)



########### Saving the results to a .csv ###########

# Define the output file path
os.makedirs("results", exist_ok=True)
output_path = os.path.join("results", "EquatorialVelocity.csv")

# Create a DataFrame to store the results
data = pd.DataFrame({
    "ROI_Center": [roi_center],
    "ROI_Size": [roi_size],
    "Temporal_Window": [temporal_window],
    "Mean_Velocity_X": [meanvx],
    "Mean_Velocity_Y": [meanvy]
})

# Check if the output file already exists
if os.path.exists(output_path):
    # Append to the existing file without writing the header again
    data.to_csv(output_path, mode='a', header=False, index=False)
else:
    # Write the DataFrame to a new CSV file
    data.to_csv(output_path, index=False)