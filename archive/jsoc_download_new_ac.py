from astropy.time import Time
import os
import datetime
from sunpy.net import Fido, attrs as a
import astropy.units as u
import numpy as np

def check_path(path, folder):
    '''
    checks if a specific directionary ('folder') in a path already exists and creates it, if not.

    Parameters
    ----------
    path:   string
            gives the path in which the wanted directionary should be in

    folder: string
            the name of the folder/directionary, that should be checked for its existence or be created, if it does not already exist
    Returns
    -------
    folder_path
        string, a combined path to the folder that was checked/created

    :Authors:
        Emily Joe Loessnitz (2024)
    '''
    folder_path = os.path.join(path, folder)
    existence = os.path.exists(folder_path)
    if not existence:
        os.makedirs(folder_path)
    return folder_path;

def check_file(path):
    if os.path.exists(path):
        os.remove(path)
        print('existing results deleted')


def downloadSDO(series, use_time_stamps, time_stamps, start_date, end_date, download_path, dt=720):

    """
    make a search query to the SDO Database to extract images of the sun in the desired timeframe.
    These are the images at --- and already correct the effects of limbdarkening.
    Checking the database might take some minutes.
    If the search-request is successful, an email will be send to your email account and the download will begin.
    The images will be downloaded as .FITS files into the directory.

    Parameters
    ----------
    t_start :  ?
            stets the start date and time in the format YYYY-MM-DDT00:00:00

    dt  :   float
            specifies the cadance between two images in hours

    days  : float?
            specifies the number of days over wich the observation shall take place

    download_path : str
            sets path and directionary for images to be downloaded to

    Returns
    -------
    numpy.ndarray
        Sorted array of data with renumbered 'id' values.

    :Authors:
        Emily Joe Loessnitz (2023)
    """
    if use_time_stamps==False:

        #steps = int(days*24 / dt)

        print(f'Preparing to download images starting from {start_date}, until {end_date}, with a {dt} second cadence.')

        print('\n\n\nChecking Database...')


        # Calculate start and end times for the search
        start_time = start_date
        end_time = end_date

        # Make the search query
        res = Fido.search(a.Time(start_time.iso, end_time.iso), a.jsoc.Series(series), a.Sample(dt*u.second) , a.jsoc.Notify(myemail))


        print(res.show('TELESCOP', 'INSTRUME', 'T_OBS'))

        #i would like to count the number of images here and double check if these equal the number of calculated steps, because sometimes there is data missing and that messes with the calculations later

        print('Does this look correct?')
        resp = input('Press [y] to download and [n] to exit\n')
        if resp == 'y':
            print('Downloading ...')
            downloaded_files = Fido.fetch(res, path=download_path+'/{file}')
        else:
            print('END')
            pass
        print('END')

    elif use_time_stamps== True:

        array_length = len(time_stamps)
        print(f'Preparing to download {array_length} images according to the provided time-stamp list...')

        # Convert string time stamps to SunPy time objects
        times = [parse_time(ts) for ts in time_stamps]

        # Create an empty list to hold individual search results
        search_results = []

        print('\n\n\nChecking Database...')

        # Perform individual searches for each time stamp
        for time in times:
            result = Fido.search(a.Time(time, time + 1 * u.minute), a.jsoc.Series(series), a.Sample(720 * u.second), a.jsoc.Notify(myemail))
            search_results.append(result)

        # Display the combined results
        for res in search_results:
            print(res)

        #i would like to count the number of images here and double check if these equal the number of calculated steps, because sometimes there is data missing and that messes with the calculations later

        print('Does this look correct?')
        resp = input('Press [y] to download and [n] to exit\n')
        if resp == 'y':
            print('Downloading ...')
            downloaded_files = Fido.fetch(*search_results, path=download_path+'/{file}')
        else:
            print('END')
            pass
        print('END')


    else:
        print('Specify if you want to download in a cadance (dt) or if you want to use an array of specific time stamps. Please set use_time_stamps to either True or False and set a cadance if False')


    # directory/folder path
    dir_path = download_path
    dirs = sorted(os.listdir(dir_path))

    # list to store files
    file_names = []

    if series== 'hmi.Ic_noLimbDark_720s':
         file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + "file_names_list.txt", "w+")
    else:
        #file_names_list = open("file_names_list.txt", "w+")
        file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + '_LimbDark_' +  "file_names_list.txt", "w+")

    # Iterate directory=
    for file_path in dirs:
        name = str(file_path)
        file_names_list.write(name + '\n')
        file_names.append(name)

    file_names_list.close()

    file_name_array = np.asarray(file_names)


    print('Download complete!')

    return FITS_path, file_name_array ;

myemail="acase2@sheffield.ac.uk"
series = 'hmi.Ic_45s'
NOAA = '0'
wavelength = 'Fe' #can be whatever, just for the folder name
mydir = ''
# not used here
time_stamps = ['2014-11-23T21:00:00', '2014-11-23T22:00:00', '2014-11-23T23:00:00']
use_time_stamps = 0

# are used here
start_date = Time('2017-09-10T15:51:00')
end_date = Time('2017-09-10T15:53:00')
dt = 45 #[seconds]

main_dir = 'NOAA_' + str(NOAA) + '-' + wavelength
new_path = check_path(mydir, main_dir)
#inside the new directionary, set up sub-dirs for images and different methods
FITS_path    = check_path(new_path, "FITS_files")


file_name = downloadSDO(series, use_time_stamps, time_stamps, start_date, end_date, FITS_path, dt)

