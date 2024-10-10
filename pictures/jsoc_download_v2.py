import datetime as dt
import os
import tarfile
import drms
import numpy as np

class client:
	_client = None

	def __init__(self, email, verbose = False):
		self._create_client(email, verbose)

	def _create_client(self, email, verbose):
		self._client = drms.Client(email = email, verbose = verbose)
		
	def create_request_string(self, series, starttime, endtime = '', wavelength = '', segment = '', period = '', cadence = ''):
		request_string = series
		if endtime: 
			if not isinstance(starttime, dt.datetime): starttime = dt.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%S')
			if not isinstance(endtime, dt.datetime): endtime = dt.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%S')
			period = abs(starttime - endtime)
		if period:
			if isinstance(period, dt.timedelta): period = str(period.total_seconds()) + 's'
			period = '/' + period
		if cadence:
			if isinstance(cadence, dt.timedelta): cadence = str(cadence.total_seconds()) + 's'
			cadence = '@' + cadence
		request_string += '[{:}{:}{:}]'.format(starttime, period, cadence)		

		if wavelength: 
			wavelength = '[{:}]'.format(wavelength)
			request_string += wavelength		
		if segment: 
			segment = '{{{:}}}'.format(segment)
			request_string += segment	
		return request_string

	def search(self, request_string, keys = ['t_obs','EXPTIME']):	
		keys = ','.join(keys)
		keys = self._client.query(request_string, key = keys)    
		return keys


	def download(self, request_string, download_dir, method = 'url-tar', protocol = 'fits', filter = None):	
		download_dir += '/'
		if not os.path.isdir(download_dir): os.makedirs(download_dir)

		if len(filter) > 0: 
			filter = np.array(filter)
			method = 'url'

		export_request = self._client.export(request_string, method = method, protocol = protocol)		

		if 'tar' in method:
			tar_downloaded = export_request.download(download_dir)['download'][0]
			with tarfile.open(tar_downloaded, "r") as tf:
				members, names = tf.getmembers(), tf.getnames()
				valid_fits = ['.fits' in filename for filename in names]
				files_downloaded = [os.path.basename(name) for name, valid in zip(names, valid_fits) if valid == True]
				members_to_extract = [member for member, valid in zip(members, valid_fits) if valid == True]	
				tf.extractall(path = download_dir, members = members_to_extract)
				os.remove(tar_downloaded)
		if not 'tar' in method:
			if np.any(filter): 
				index = np.where(filter == True)[0] 
				files_downloaded = export_request.download(download_dir, index = index)['download'][index]
			else: 
				files_downloaded = export_request.download(download_dir)['download']
			files_downloaded = [os.path.basename(file) for file in files_downloaded]
			for i, file in enumerate(files_downloaded):
				name, ending = os.path.splitext(file)
				if ending != '.fits':
					os.replace(download_dir + file, download_dir + name)
					files_downloaded[i] = name			
		return files_downloaded






