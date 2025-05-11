import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Find delta_t for which d = 1.25 Mm
target_distance = 1.25  # Target mean distance in Mm

# Define the path to the output directory
output_dir = r'C:/Users/AlexC/source/repos/solar_granulations/output'

# Initialise a dictionary to store data
file_data = {}

# Define unique markers for each file
markers = ['o', 'x', 'D']  # Add more markers if needed
marker_index = 0

# Define custom labels for the files
custom_labels = {
    '3.5 minute threshold.txt': '5 minutes',
    '1.135s threshold.txt': '135 seconds',
    '2. No threshold.txt': 'No thresholding'
}


# Loop through the .txt files in the output directory
for file_name in os.listdir(output_dir):

    file_path = os.path.join(output_dir, file_name)
    
    # Extract the data, skip the first line
    data = pd.read_csv(file_path, sep='\t', skiprows=1, names=['Temporal_Window', 'Mean_Distance'])

    # Convert Mean_Distance from pixels to Mm (1 pixel = 0.725 Mm)
    data['Mean_Distance'] = data['Mean_Distance'] * 0.725
    
    # Store the data in the dictionary with the file name as the key
    file_data[file_name] = data

    # Interpolate the data
    interpolation = interp1d(data['Mean_Distance'], data['Temporal_Window'], kind='linear', fill_value="extrapolate")
    
    # Calculate delta_t for the target distance
    delta_t = interpolation(target_distance)
    print(f"For file '{file_name}', delta_t for d = {target_distance} Mm is approximately {delta_t:.2f} seconds.")



# Plot the data from each file as a separate line
plt.figure(figsize=(10, 10))
for file_name, data in file_data.items():

    # Use a unique marker for each file
    marker = markers[marker_index % len(markers)]
    marker_index += 1

    # Use custom labels if available, otherwise use the file name
    label = custom_labels.get(file_name, file_name)

    plt.plot(
        data['Temporal_Window'],
        data['Mean_Distance'],
        marker=marker,
        label=label,
        color='black',
        linestyle='--',  # Solid line for all files
    )

plt.xscale('log', base=10)

# Set x-axis limits to start at 10^2 and end at 2*10^3
plt.xlim(10**2, 2 * 10**3)

# Explicitly set x-axis ticks
plt.xticks([10**2, 10**3, 2 * 10**3], labels=["$10^2$", "$10^3$", "$2 \\times 10^3$"])

plt.title('Mean Distance between Measured Velocities vs Temporal Window')
plt.xlabel('Temporal Window (s)')
plt.ylabel('Mean Distance (Mm)')
plt.legend(title='Thresholds')
plt.grid()

os.makedirs('plots', exist_ok=True)
output_path = os.path.join('plots', 'mean_distance_vs_temporal_window.png')
plt.savefig(output_path)
print(f'Mean Distance plot saved to {output_path}')