import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend for matplotlib

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

# Define the file folder where the FITS files are stored
file_folder = 'test/testdata'
file_list = sorted(os.listdir(file_folder))  # List the files in the folder
file_index = 0

# Storing FITS image as SunPy Map object
sunpy_map = sm.Map(f'{file_folder}/{file_list[file_index]}')  # Create a SunPy map object from the first image
sunpy_map_rotated = sunpy_map.rotate(order = 3)  # Rotate the SunPy map object using bi-cubic iterpolation so the solar north is at the top of the image

# Display the SunPy map using Matplotlib and save the figure
plt.figure(figsize=(10, 10))
sunpy_map_rotated.plot()
plt.title("SunPy Map, 2024/08/15, 08:00:45")
output_path = os.path.join('test', 'plots', 'sunpy_map.png')
plt.savefig(output_path, dpi=1000)