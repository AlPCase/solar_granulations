### Load env, change to folder, start python
# conda activate your_work_enviroment
# cd /home/your_username/filepath/
# jsoc_download_v2.py should be in the same directory
# need dataetime and os modules for python
# python

### 
import jsoc_download_v2 as jd
import datetime as dt
import os

mystarttime='2017-09-10T15:51:00' #**************************************
myendtime='2017-09-10T15:52:00'   #**************************************
mystart_time = dt.datetime.strptime(mystarttime,'%Y-%m-%dT%H:%M:%S')
myend_time = dt.datetime.strptime(myendtime, '%Y-%m-%dT%H:%M:%S')
myseries='hmi.Ic_45s'             #**************************************
mywavelengths="6173.3"               #**************************************  
myemail="acase2@sheffield.ac.uk" #**************************************


# Create a client instance
my_client = jd.client(myemail, verbose=True)

# Create a request string
request_string = my_client.create_request_string(myseries, mystart_time, myend_time, wavelength=mywavelengths, segment="image")

# Search using the request string
download_dir = 'C:/Users/AlexC/source/repos/solar_granulations/imagedata'  #*******************************************
search_result = my_client.search(request_string)
print(search_result)

# Download using the request string
downloaded_files = my_client.download(request_string, download_dir, filter = search_result['EXPTIME'] < 1)
print(downloaded_files)

