import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from sunpy.net import Fido, attrs as a
import sunpy.map
from astropy.io import fits
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
#file_path = os.path.join(file_folder, "xyVelocities.csv")

# Load the data from the CSV file
data = pd.read_csv(file_folder)

# Extract y velocity and roi_center[1]
y_velocity = data['Mean_Velocity_Y']
x_velocity = data['Mean_Velocity_X']
roi_center = data['ROI_Center'].apply(ast.literal_eval)
roi_center_y = [coord[1] for coord in roi_center]

# Convert using FITS header metadata
obstime = hmi_map.date.to_datetime().strftime("%Y-%m-%dT%H:%M:%S")  # Get time from FITS
solar_latitudes = [
    arcsec_to_solar_latitude(y, observer, obstime) 
    for y in roi_center_y
]

print(f"Y Velocities: {y_velocity}")
print(f"ROI Centers: {roi_center}")
print(f"ROI Center Y: {roi_center_y}")

# Filter data for solar latitudes within ±60 degrees
filtered_latitudes = []
filtered_velocities_y = []
filtered_velocities_x = []

for lat, vely, velx in zip(solar_latitudes, y_velocity, x_velocity):
    if abs(lat) <= 60:
        filtered_latitudes.append(lat)
        filtered_velocities_y.append(vely)
        filtered_velocities_x.append(velx)


# Update the data for plotting
solar_latitudes = filtered_latitudes
y_velocity = filtered_velocities_y
x_velocity = filtered_velocities_x
# Plot the data
plt.figure(figsize=(8, 6))
plt.plot(solar_latitudes, y_velocity, color='black')
plt.ylim(0, 1.3)
plt.xlabel('Heliographic Latitude (degrees)')
plt.ylabel('Y Velocity (km/s)')
plt.title('Differential Rotation\n Y Velocity vs Heliographic Latitude')
plt.legend()
plt.grid(True)

# Define the output directory
os.makedirs('plots', exist_ok=True)
output_file = os.path.join('plots', 'y_velocity_vs_roi_center_y.png')
plt.savefig(output_file)
print(f'Plot saved to {output_file}')