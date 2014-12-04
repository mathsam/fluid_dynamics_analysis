from scipy.io import netcdf
filename = '/archive/Junyi.Chai/QG_exp/Aug13/Sc=1.6/0.1_kd160_resox2/HISTORY/history007.nc'

f = netcdf.netcdf_file(filename,'r')
vor_p     = f.variables['PV_anomaly']
vor_last  = vor_p[-1]
vor_upper = vor_last[0,:,:]
vor_lower = vor_last[1,:,:]

ucomp_p   = f.variables['ug']
ucomp     = ucomp_p[-1]
## plotting
import numpy as np
import matplotlib.pyplot as plt

PI = 3.1415926
N  = 512  # 512
x = np.linspace(-PI,PI,N)
y = np.linspace(-PI,PI,N)
#----- draw vorticity -------------
fig = plt.figure()
ax1 = fig.add_subplot(221)
im1 = ax1.pcolor(x,y,vor_upper,cmap=plt.cm.bone)
#plt.axis('equal')
plt.title('upper')
ax1.set_aspect('equal')
ax1.set_xlim([-PI,PI])
ax1.set_ylim([-PI,PI])
fig.colorbar(im1)


#plt.axes(axes[0][0])
#c_range = 2*vor_upper.std()
#plt.clim([-c_range,c_range])

ax2 = fig.add_subplot(222)
im2 = ax2.pcolor(x,y,vor_lower,cmap=plt.cm.bone)


#plt.clim([-c_range,c_range])
#axarr[1].colorbar()
plt.title('lower')
ax2.set_aspect('equal')
ax2.set_xlim([-PI,PI])
ax2.set_ylim([-PI,PI])
fig.colorbar(im2)

#-------- draw velocity ------------------
ax1 = fig.add_subplot(223)
im1 = ax1.pcolor(x,y,ucomp[0,:,:],cmap=plt.cm.bone)
#plt.axis('equal')
ax1.set_aspect('equal')
ax1.set_xlim([-PI,PI])
ax1.set_ylim([-PI,PI])
fig.colorbar(im1)


#plt.axes(axes[0][0])
#c_range = 2*vor_upper.std()
#plt.clim([-c_range,c_range])

ax2 = fig.add_subplot(224)
im2 = ax2.pcolor(x,y,ucomp[1,:,:],cmap=plt.cm.bone)


#plt.clim([-c_range,c_range])
#axarr[1].colorbar()
ax2.set_aspect('equal')
ax2.set_xlim([-PI,PI])
ax2.set_ylim([-PI,PI])
fig.colorbar(im2)

plt.show()