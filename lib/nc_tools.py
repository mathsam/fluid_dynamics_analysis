from scipy.io import netcdf
import os
import re
import numpy as np
import logging
from list_tools import ComparableList
from mmap_array import MmapArray


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
_log_handler = logging.StreamHandler()
_formater = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
_log_handler.setFormatter(_formater)
logger.addHandler(logging.StreamHandler())

def ncread(filedir, filename_regexp, var_name):
    """
    A wrapper for NetCDFChain
    If there is only one file, return netcdf_varialbe;
    otherwise, return NetCDFChain object.
    """
    files  = os.listdir(filedir)
    regexp = re.compile(filename_regexp)
    filtered_files = []
    for each_file in files:
        if (regexp.search(each_file)):
            filtered_files.append(each_file)

    if (len(filtered_files) == 1):
        if filedir[-1] != '/':
            filedir += '/'
        filename = filedir + filtered_files[0]
        fh = netcdf.netcdf_file(filename,'r',mmap=True)
        return np.array(fh.variables[var_name][:])
    else:
        return NetCDFChain(filedir, filename_regexp, var_name)
            

class NetCDFChain(object):
    """Manage multiple NetCDF files from restart runs as a single
    file transparently. The syntax is the same as numpy array. 

    Sample usage:
        filedir  = '/archive/Junyi.Chai/QG_exp/Nov28_kfe-2_qg'
        filename = r'Nov28_kfe-2_qg_energy_seg[0-9]+'
        f        = nc_chains(filedir,filename,'ke')
        ke       = f[:]  # combine all times together
        ke1      = f[0]  # get first element (time_dim = 0 for these files)
        ke2      = f[-1] # get the last element
        ke3      = f[1:2,0:5,:] # if the variable are multi-dimensional

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
    def __new__(cls, filedir, filename_regexp, var_name, time_dim = 0, 
                     last_n_files = None):
        """
        If there is only one file, return netcdf_varialbe;
        otherwise, return NetCDFChain object.
        """
        files  = os.listdir(filedir)
        regexp = re.compile(filename_regexp)
        filtered_files = []
        for each_file in files:
            if (regexp.search(each_file)):
                filtered_files.append(each_file)
    
        if (len(filtered_files) == 1):
            if filedir[-1] != '/':
                filedir += '/'
            filename = filedir + filtered_files[0]
            logger.info("Only one file matches: %s" %filename)
            fh = netcdf.netcdf_file(filename,'r',mmap=True)
            return MmapArray(fh.variables[var_name][:], fh)
        else:
            return object.__new__(cls)     

    def __init__(self, filedir, filename_regexp, var_name, time_dim = 0,
                       last_n_files = None):
        """filedir: string
                directory where the files are stored.
        filename_regexp: string, raw string
                 regular expression to match for filenames.
                 If just want read one file, simply use filename.
        var_name: string
              variable name to read in
        time_dim: {0, 1, 2, ...}, optional
              the dimension that one wants to combine multiple files into
              a single file.
        last_n_files: select the last n files from those that match the regexp      
        """

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

        if (len(filtered_files) <= 1):
            raise ValueError('less than two files found')

        # get numbers in filenames and sort according to it
        regexp_num = re.compile(r'[0-9]+')
        nums_in_filtered_files = []
        for each_file in filtered_files:
                filename = filedir + each_file
                #add file only if the target variable exists in it
                #however, this is very slow on tape system!
                #try:
                #    if(netcdf.netcdf_file(filename,'r').variables[var_name]):
                #       nums_in_filtered_files.append(regexp_num.findall(each_file))
                #except KeyError:
                #    logger.info("Warning: %s does not contain %s\n", 
                #                each_file, var_name)
                nums_in_filtered_files.append(regexp_num.findall(each_file))

        self.sorted_files = [each_file for (each_file, tmp) in
          sorted(zip(filtered_files, nums_in_filtered_files),
             key = lambda li : ComparableList(li[1]))]
             
        if last_n_files:
            self.sorted_files = self.sorted_files[-last_n_files:]

        self.time_steps_each_file = []
        self.time_steps_before    = [0,]
        self.ndim = None
        vshape     = None
        axes       = None
        for each_file in self.sorted_files:
            filename = filedir + each_file
            var = netcdf.netcdf_file(filename,'r').variables[var_name]
            time_len = var.shape[time_dim]
            self.time_steps_each_file.append(time_len)
            total_time_steps = self.time_steps_before[-1] + time_len
            self.time_steps_before.append(total_time_steps)
            if(self.ndim):
                if(self.ndim != len(var.shape)):
                    raise Exception("change in variable dimensions")
            else: 
                self.ndim = len(var.shape)
            if(vshape):
                # num of time steps can change; other shape cannot change
                if(var.shape[:time_dim] + var.shape[time_dim+1:] != vshape):
                    raise Exception("change in variable shape")
            else:
                vshape = var.shape[:time_dim] + var.shape[time_dim+1:]
            if(axes):
                if(var.dimensions != axes):
                    raise Exception("change in variable dimensions (axes)")
            else:
                axes = var.dimensions

        self.total_time_steps = total_time_steps
        self.shape = vshape[:time_dim] + (total_time_steps,) + vshape[time_dim:]
        self.dimensions = axes

        logger.info("Files contained:\n"
                  + self.sorted_files[0] + '\n...\n' + self.sorted_files[-1]
                  + "\n%d number of time steps in total\n",total_time_steps)
                  
    def __repr__(self):
        return "NetCDFChain object\nContains\n%s\n...\n%s\n%d time steps in total\n%s as axes (dimensions)" %(self.sorted_files[0], self.sorted_files[-1], 
                 self.total_time_steps,
                 self.dimensions)


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
        
        if index.start is not None and index.start < 0:
            if index.start >= -self.total_time_steps:
                index = slice(index.start + self.total_time_steps, 
                    index.stop, index.step)
            else:
                raise IndexError("%d is out of range" %index.start)
                
        if index.stop is not None and index.stop < 0:
            if index.stop >= -self.total_time_steps:
                index = slice(index.start, 
                    index.stop + self.total_time_steps, index.step)
            else:
                raise IndexError("%d is out of range" %index.stop)
            

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

        logger.info("Index [%s:%s:%s] is out of range\n",str(index.start), 
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
            f   = netcdf.netcdf_file(filename,'r',mmap=True)
            var = f.variables[self.var_name]
            #var = netcdf.netcdf_file(filename,'r').variables[self.var_name]
            if(combined_var is None):
                combined_var = np.array(var[each_slice])
            else:
                tmp_var      = np.array(var[each_slice])
                combined_var = np.concatenate((combined_var,tmp_var),
                                              axis = self.time_dim)
        return combined_var


    def __getitem__(self,index):
        if isinstance(index, int):
            if (self.time_dim == 0):
                file_index, time_slice = self._time_to_file(index)
                logger.info("file = %s, local time index = %d\n",
                            self.sorted_files[file_index], time_slice)
                filename = self.filedir + self.sorted_files[file_index]
                fh = netcdf.netcdf_file(filename,'r',mmap=True)
                var = fh.variables[self.var_name]
                return np.array(var[time_slice])
            else:
                one_slice = [slice(None)]*self.ndim
                one_slice[0] = index
                slice_list   = [tuple(one_slice)] * len(self.sorted_files)
                return self._combine_slices(slice_list)

        if isinstance(index, slice):
            if(self.time_dim != 0):
                one_slice = [slice(None)]*self.ndim
                one_slice[0] = index
                slice_list   = [tuple(one_slice)] * len(self.sorted_files)
                return self._combine_slices(slice_list)
            else:
                slice_map, file_id_list = self._timeslice_to_file(index)
                NULL_slice = [slice(None)]*(self.ndim-1)
                slice_list = [tuple([each_slice] + NULL_slice) for each_slice 
                              in slice_map]
                return self._combine_slices(slice_list,file_id_list)
                
        if index is Ellipsis:
            index = slice(None)
            return self[index]

        if isinstance(index,tuple):
            if len(index) > self.ndim:
                raise KeyError('index has more dimensions than the array')
            if Ellipsis in index:
                first_ellip = index.index(Ellipsis)
                index = index[0:first_ellip] + (slice(None),)*(self.ndim 
                    - len(index) + 1) + index[first_ellip+1:]
                index = [slice(None) if i is Ellipsis else i for i in index]
                return self[tuple(index)]

            one_slice = []
            for each_index in index:
                if isinstance(each_index, int):
                    one_slice += [each_index]
                elif isinstance(each_index, slice):
                    one_slice += [each_index]
            one_slice += [slice(None)]*(self.ndim - len(index))
    
            # if time dim is outside index, need all time steps
            if(self.time_dim > len(index)-1):
                slice_list = [tuple(one_slice)] * len(self.sorted_files)
                return self._combine_slices(slice_list)
    
            time_index = index[self.time_dim]
            if isinstance(time_index,int):
                file_index, time_slice = self._time_to_file(time_index)
                filename = self.filedir + self.sorted_files[file_index]
                fh = netcdf.netcdf_file(filename,'r',mmap=True)
                var = fh.variables[self.var_name]
                one_slice[self.time_dim] = time_slice
                return np.array(var[tuple(one_slice)])
            elif isinstance(time_index,slice):
                slice_map, file_id_list = self._timeslice_to_file(time_index)
                slice_list = []
                for each_slice in slice_map:
                    one_slice[self.time_dim] = each_slice
                    slice_list.append(tuple(one_slice))
                return self._combine_slices(slice_list,file_id_list)
            else:
                raise Exception("Unknown situation")

class CreateRestartNC(object):
    def __init__(self, nc_name='restart.nc', kx=1023, ky=512, z=1, ztracer=1,version=1):
        """Initialize a NetCDF file with axis kx, ky, z, ztracer of specified
        length
        This intends to be the restart file so time dimension only has length 1
        
        Inputs:
            nc_name: filename for the NetCDF file
            kx, ky, z, ztracer: axis length for those axes
                Not that if input is None, does not create this axis
            version: version of netcdf to read / write, where 1 means Classic 
                format and 2 means 64-bit offset format. The QG model uses classic
        """
        self.nc_name = nc_name
        self.f = netcdf.netcdf_file(nc_name, 'w', version=version)
        self.kx = kx
        self.ky = ky
        self.z  = z
        self.f.createDimension('time_step', 1)
        self.f.createDimension('real_and_imag', 2)
        self.f.createDimension('kx', kx)
        self.f.createDimension('ky', ky)
        self.f.createDimension('z',  z)
        if ztracer:
            self.ztracer = ztracer
            self.f.createDimension('ztracer', ztracer)
    
    def write_var(self, var_name, var, is_tracer=False):
        """write a variable `var` with name `var_name` into the file
        
        Inputs:
            var_name: string
            var: complex numpy array with shape (t=1, ky, kx, z|ztracer)
        """
        if var.ndim != 4:
            raise ValueError('dimension mismatch with requirement')
        if var.shape != (1, self.ky, self.kx, self.z) or var.shape != (
                         1, self.ky, self.kx, self.ztracer):
            raise ValueError('shape mismatch with initialization')
        if not is_tracer:
            empty_var = self.f.createVariable(var_name, 'float', 
                ('time_step','real_and_imag','ky','kx','z'))
        else:
            empty_var = self.f.createVariable(var_name, 'float', 
                ('time_step','real_and_imag','ky','kx','ztracer'))
        empty_var[:,0,:,:,:] = var.real
        empty_var[:,1,:,:,:] = var.imag
    
    def __del__(self):
        self.f.flush()
        self.f.close()


if __name__ == "__main__":
    nc_name = 'test_restart.nc'
    psi = np.random.normal(0., 1., (1, 512, 1023, 1)) \
          + 1j*np.random.normal(0., 1., (1, 512, 1023, 1))
    f = CreateRestartNC(nc_name)
    f.write_var('psi',psi)
    del f
