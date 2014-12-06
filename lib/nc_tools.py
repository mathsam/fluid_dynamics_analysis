from scipy.io import netcdf
import os
import re
import numpy as np
from list_tools import ComparableList

class nc_chains(object):
    """Manage multiple NetCDF files from restart runs as a single
    file transparently. 

    Sample usage:
        filedir  = '/archive/Junyi.Chai/QG_exp/Nov28_kfe-2_qg'
        filename = r'Nov28_kfe-2_qg_energy_seg[0-9]+'
        f        = nc_chains(filedir,filename,'ke')
        ke       = f[:]  # combine all times together
        ke1      = f[0]  # get time zero (time_dim = 0 for these files)
        ke2      = f[-1] # get the last element

    The files are automatically sequenced by the numbers
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

        self.var_name = var_name
        self.time_dim = time_dim
        if filedir[-1] != '/':
            filedir += '/'
            self.filedir = filedir

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
                filename = filedir + each_file
                #add file only if the target variable exists in it
                try:
                    if(netcdf.netcdf_file(filename,'r').variables[var_name]):
                       nums_in_filtered_files.append(regexp_num.findall(each_file))
                except KeyError:
                    print "Warning: %s does not contain %s" %(each_file, var_name)

        self.sorted_files = [each_file for (each_file, tmp) in
          sorted(zip(filtered_files, nums_in_filtered_files),
             key = lambda li : ComparableList(li[1]))]

        self.time_steps_each_file = []
        self.time_steps_before    = [0,]
        self.num_dims = None
        for each_file in self.sorted_files:
            filename = filedir + each_file
            var = netcdf.netcdf_file(filename,'r').variables[var_name]
            time_len = var.shape[time_dim]
            self.time_steps_each_file.append(time_len)
            total_time_steps = self.time_steps_before[-1] + time_len
            self.time_steps_before.append(total_time_steps)
            if(self.num_dims):
                if(self.num_dims != len(var.shape)):
                    raise Expection("change in variable dimensions")
            else: self.num_dims = len(var.shape)

        self.total_time_steps = total_time_steps

        print "Files contained:"
        print self.sorted_files[0] + '\n...\n' + self.sorted_files[-1]
        print "%d number of time steps in total" %total_time_steps


    def _time_to_file(self,time):
        """map a time (which is an int) to the corresponding 
        file and the time slice inside this file. Indexing similar to that of a 
        list. For example,
            0 means the first time
            -1 means the last time"""
        if(time >= 0):
            if(time >= self.total_time_steps):
                raise Exception("Time indexing is out of range")
            file_index = next(x[0] for x in enumerate(self.time_steps_before[1:])
                              if time < x[1])
            time_slice = time - self.time_steps_before[file_index]
            return file_index, time_slice
        elif(time >= -self.total_time_steps):
            return self._time_to_file(time + self.total_time_steps)
        else:
            raise Exception("Time indexing is out of range")


    def _timeslice_to_file(self, index):
        """Given a time range, represented by index (which is a slice obj), 
        return the corresponding slice in each file and the file ids
        
        index: slice object
               For example, slice(0, 10, 1)"""

        slice_list   = []
        file_id_list = []
        if(not index.start):
            index = slice(0, index.stop, index.step)

        for i in (x[0] for x in enumerate(self.sorted_files)):
            # the range of indexes in file i is 
            # [self.time_steps_before[i],
            #  self.time_steps_before[i] + self.time_steps_each_file[i] -1]
            # where [] is a closed range

            if(index.start <= self.time_steps_before[i]
                            + self.time_steps_each_file[i] -1):
                one_slice = slice(None,None,index.step)
                if(index.start >= self.time_steps_before[i]):
                    one_slice = slice(index.start - self.time_steps_before[i],
                                            one_slice.stop, one_slice.step)
                else:
                    one_slice = slice(0, one_slice.stop, one_slice.step)

                if(index.stop is not None and 
                   index.stop <= self.time_steps_before[i] 
                               + self.time_steps_each_file[i] -1):
                    one_slice = slice(one_slice.start,
                                      index.stop - self.time_steps_before[i],
                                      one_slice.step)
                    slice_list.append(one_slice)
                    file_id_list.append(i)
                    return slice_list, file_id_list
                else: # include the case index.stop == None
                    one_slice = slice(one_slice.start, None, one_slice.step)
                    slice_list.append(one_slice)
                    file_id_list.append(i)
        if slice_list:
            return slice_list, file_id_list

        print "Index [%s:%s:%s] is out of range" %(str(index.start), 
                                                   str(index.stop), 
                                                   str(index.step))
        raise Exception("Indexing returns no files")


    def _combine_slices(self,slice_list,file_id_list = None):
        """given a list of slices, for example[(0:2,0:5), (0:2,0:5)], and a list
        of the indexes in self.sorted_files, then combine the variable in with 
        these slices/indexes from the list_lists together

        slice_list  : a list of tuples
        file_id_list: a list of integers | None
                      if None, means all the files"""
        combined_var = None
        if(file_id_list is None):
            file_id_list = [i for i,tmp in enumerate(self.sorted_files)]

        for each_file_id, each_slice in zip(file_id_list, slice_list):
            filename = self.filedir + self.sorted_files[each_file_id]
            f   = netcdf.netcdf_file(filename,'r',mmap=False)
            var = f.variables[self.var_name]
            #var = netcdf.netcdf_file(filename,'r').variables[self.var_name]
            if(combined_var is None):
                combined_var = var[each_slice]
            else:
                tmp_var      = var[each_slice]
                combined_var = np.concatenate((combined_var,tmp_var),
                                              axis = self.time_dim)
        return combined_var


    def __getitem__(self,index):
        if isinstance(index, int):
            if (self.time_dim == 0):
                file_index, time_slice = self._time_to_file(index)
                print "file = %s, local time index = %d" \
                       %(self.sorted_files[file_index], time_slice)
                filename = self.filedir + self.sorted_files[file_index]
                var = netcdf.netcdf_file(filename,'r',mmap=False).variables[self.var_name]
                return var[time_slice] 
            else:
                one_slice = [slice(None)]*self.num_dims
                one_slice[0] = index
                slice_list   = [tuple(one_slice)] * len(self.sorted_files)
                return self._combine_slices(slice_list)

        if isinstance(index, slice):
            if(self.time_dim != 0):
                one_slice = [slice(None)]*self.num_dims
                one_slice[0] = index
                slice_list   = [tuple(one_slice)] * len(self.sorted_files)
                return self._combine_slices(slice_list)
            else:
                slice_map, file_id_list = self._timeslice_to_file(index)
                NULL_slice = [slice(None)]*(self.num_dims-1)
                slice_list = [tuple([each_slice] + NULL_slice) for each_slice 
                              in slice_map]
                return self._combine_slices(slice_list,file_id_list)

        if isinstance(index,tuple):
           one_slice = []
           for each_index in index:
               if isinstance(each_index, int):
                   one_slice += [each_index]
               elif isinstance(each_index, slice):
                   one_slice += [each_index]
           one_slice += [slice(None)]*(self.num_dims - len(index))

            # if time dim is outside index, need all time steps
           if(self.time_dim > len(index)-1):
               slice_list = [tuple(one_slice)] * len(self.sorted_files)
               return self._combine_slices(slice_list)

           time_index = index[self.time_dim]
           if isinstance(time_index,int):
               file_index, time_slice = self._time_to_file(time_index)
               filename = self.filedir + self.sorted_files[file_index]
               var = netcdf.netcdf_file(filename,'r',mmap=False).variables[self.var_name]
               one_slice[self.time_dim] = time_slice
               return var[tuple(one_slice)] 
           elif isinstance(time_index,slice):
               slice_map, file_id_list = self._timeslice_to_file(time_index)
               slice_list = []
               for each_slice in slice_map:
                   one_slice[self.time_dim] = each_slice
                   slice_list.append(tuple(one_slice))
               return self._combine_slices(slice_list,file_id_list)
           else:
               raise Exception("Unknown situation")


if __name__ == "__main__":
    filename_regexp = r'Nov28_qg_seg[0-9]+'
    filedir = '/archive/Junyi.Chai/QG_exp/Nov28_qg'
    f = nc_chains(filedir,filename_regexp,'psi')
    a = f[49:51]
    print a.shape
