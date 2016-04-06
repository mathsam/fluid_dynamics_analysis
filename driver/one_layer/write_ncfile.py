x = np.zeros((1024, 1024))
y = np.zeros((1024, 1024))
z = np.zeros((1024, 1024))
for i, phi in enumerate(lons/180*np.pi+np.pi):
    print i
    for j, theta in enumerate(lats/90*np.pi+np.pi):
        x[j,i] = (4 + 2*np.cos(theta))*np.cos(phi)
        y[j,i] = (4 + 2*np.cos(theta))*np.sin(phi)
        z[j,i] = 2*np.sin(theta)
##
from netCDF4 import Dataset
import numpy as np

root_grp = Dataset('torus.nc', 'w', format='NETCDF4')
root_grp.description = 'vorticity data'

# dimensions
root_grp.createDimension('lat', 1024)
root_grp.createDimension('lon', 1024)

# variables
latitudes = root_grp.createVariable('lat', 'f4', ('lat',))
longitudes = root_grp.createVariable('lon', 'f4', ('lon',))
vor = root_grp.createVariable('vor', 'f4', ('lat', 'lon',))
torusx = root_grp.createVariable('torusx', 'f4', ('lat', 'lon',))
torusy = root_grp.createVariable('torusy', 'f4', ('lat', 'lon',))
torusz = root_grp.createVariable('torusz', 'f4', ('lat', 'lon',))

# data
lats =  np.linspace(-90, 90, 1024)
lons =  np.linspace(-180, 180, 1024)
latitudes[:] = lats
longitudes[:] = lons
vor[:,:] = vorg2d[:,:]

torusx[:,:] = x[:,:]
torusy[:,:] = y[:,:]
torusz[:,:] = z[:,:]        

root_grp.close()