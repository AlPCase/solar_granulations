import os
import astropy.io.fits as fits



file_folder = "NOAA_0-Fe/FITS_files"  # Define the file folder where the FITS files are stored ############# Should make this a global variable
file_list = sorted(os.listdir(file_folder))  # List the files in the folder


hdr = fits.open(f'{file_folder}/{file_list[1]}')[1].header
print(hdr)
