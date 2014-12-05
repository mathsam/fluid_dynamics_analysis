from scipy.io import netcdf
import os
import re
from list_tools import ComparableList

class nc_chains(object):
    """Manage multiple NetCDF files from restart runs as a single
    file transparently. The files are automatically sequenced by the numbers
    in the file names. If multiple numbers appear in filename, the files are
    sorted by the first appearance of number first, and then the second, and
    so on.

    For example:
       Nov27_seg1.nc
       Nov27_seg2.nc
       Nov28_seg0.nc
       Nov28_seg1.nc
       ...
    """
    def __init__(self, filedir, filename_regexp, var_name, time_dim = 0):
	"""filedir: string
		    directory where the files are stored.
	filename_regexp: string, raw string
			 regular expression to match for filenames.
			 If just want read one file, simply use filename.
	var_name: string
		  variable name to read in
	time_dim: {0, 1, 2, ...}, optional
		  the dimension that one wants to combine multiple files into
		  a single file."""
	if filedir[-1] != '/':
		filedir += '/'

	files  = os.listdir(filedir)
	regexp = re.compile(filename_regexp)
	filtered_files = []
	for each_file in files:
	    if (regexp.search(each_file)):
		filtered_files.append(each_file)

	if (len(filtered_files) == 1):
	    filename = filedir + filtered_files[0]
	    return netcdf.netcdf_file(filename,'r').variables[var_name]

        # get numbers in filenames and sort according to it
	regexp_num = re.compile(r'[0-9]+')
	nums_in_filtered_files = []
	for each_file in filtered_files:
	    nums_in_filtered_files.append(regexp_num.findall(each_file))

	self.sorted_files = [each_file for (each_file, tmp) in
	  sorted(zip(filtered_files, nums_in_filtered_files),
		 key = lambda li : ComparableList(li[1]))]
        print "Files contained:"
        print self.sorted_files[0] + '\n...\n' + self.sorted_files[-1]

        self.time_steps_each_file = []
        self.time_steps_before    = [0,]
        for each_file in self.sorted_files:
            filename = filedir + each_file
            var = netcdf.netcdf_file(filename,'r').variables[var_name]
            time_len = var.shape[time_dim]
            self.time_steps_each_file.append(time_len)
            total_time_steps = self.time_steps_before[-1] + time_len
            self.time_steps_before.append(total_time_steps) 
