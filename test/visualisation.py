import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

# Define the path to the output directory
output_dir = r'C:/Users/AlexC/source/repos/solar_granulations/output'

# Initialise a dictionary to store data
file_data = {}

# Loop through the .txt files in the output directory
for file_name in os.listdir(output_dir):

    file_path = os.path.join(output_dir, file_name)
    
    # Extract the data, skip the first line
    data = pd.read_csv(file_path, sep='\t', skiprows=1, names=['Temporal_Window', 'Mean_Distance'])
    
    # Store the data in the dictionary with the file name as the key
    file_data[file_name] = data

# Plot the data from each file as a separate line
plt.figure(figsize=(10, 10))
for file_name, data in file_data.items():
    plt.plot(
        data["Temporal_Window"],
        data["Mean_Distance"],
        marker='o',
        label=file_name
    )

#plt.xscale('log', base=10)

plt.title("Mean Distance between Measured Velocities vs Temporal Window")
plt.xlabel("Temporal Window (s)")
plt.ylabel("Mean Distance (pixels)")
plt.legend(title="Files")
plt.grid()

os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'mean_distance_vs_temporal_window.png')
plt.savefig(output_path)
print(f"Mean Distance plot saved to {output_path}")