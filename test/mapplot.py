import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
import sunpy.map
import ast
import os

# --------------------------------------------------
# Function to convert arcsec to solar latitude
# --------------------------------------------------
def arcsec_to_solar_latitude(y_arcsec, observer, obstime):
    """
    Convert y-coordinate (arcseconds) to solar latitude (degrees).
    
    Args:
        y_arcsec (float): Helioprojective y-coordinate in arcseconds.
        observer (SkyCoord): Observer location (e.g., SDO).
        obstime (str): Observation time (from FITS header).
    
    Returns:
        float: Solar latitude in degrees.
    """
    coord = SkyCoord(
        Tx=0*u.arcsec,  # Assuming x=0 for latitude conversion
        Ty=y_arcsec*u.arcsec,
        frame=frames.Helioprojective,
        obstime=obstime,
        observer=observer
    )
    heliographic = coord.transform_to(frames.HeliographicStonyhurst)
    return heliographic.lat.to(u.deg).value


# Load your FITS file to get metadata
fits_path = r'C:\Users\AlexC\source\repos\solar_granulations\NOAA_0-Fe\FITS_files\hmi.ic_45s.20170910_155145_TAI.2.continuum.fits'  # Replace with actual FITS file
hmi_map = sunpy.map.Map(fits_path)
hmi_map = hmi_map.rotate(order=3)

# Get observer metadata from FITS header
observer = hmi_map.observer_coordinate

# Define path to the CSV file
file_folder = r"C:\Users\AlexC\source\repos\solar_granulations\test\results\xyVelocities.csv"

# Load the data from the CSV file
data = pd.read_csv(file_folder)

# Extract y velocity and roi_center[1]
y_velocity = data['Mean_Velocity_Y']
roi_center = data['ROI_Center'].apply(ast.literal_eval)
roi_center_y = [coord[1] for coord in roi_center]

# Convert using FITS header metadata
obstime = hmi_map.date.to_datetime().strftime("%Y-%m-%dT%H:%M:%S")  # Get time from FITS
solar_latitudes = [
    arcsec_to_solar_latitude(y, observer, obstime) 
    for y in roi_center_y
]

# Filter data for solar latitudes within ±60 degrees
filtered_latitudes = []
filtered_velocities_y = []

for lat, vely in zip(solar_latitudes, y_velocity):
    if abs(lat) <= 60:
        filtered_latitudes.append(lat)
        filtered_velocities_y.append(vely)

# Update the data for plotting
solar_latitudes = filtered_latitudes
y_velocity = filtered_velocities_y

# Create a solar disk plot
plt.figure(figsize=(8, 8))
ax = plt.gca()

# Draw the solar disk
circle = plt.Circle((0, 0), 1, color='orange', alpha=0.5, label='Solar Disk')
ax.add_artist(circle)

# Normalize the latitudes to fit on the solar disk
normalized_latitudes = [lat / 90 for lat in solar_latitudes]  # Scale latitudes to [-1, 1]

# Plot the Y velocities as arrows
arrow_scale = 0.5  # Scale factor for arrow length
for lat, vel in zip(normalized_latitudes, y_velocity):
    plt.quiver(0, lat, vel * arrow_scale, 0, angles='xy', scale_units='xy', scale=1, color='black')

# Set plot limits and labels
plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.title('Y Velocities Overlaid on Solar Disk')
plt.gca().set_aspect('equal', adjustable='box')

# Save the plot
os.makedirs('plots', exist_ok=True)
output_file = os.path.join('plots', 'y_velocity_on_solar_disk.png')
plt.savefig(output_file)
print(f'Plot saved to {output_file}')