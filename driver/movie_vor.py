from scipy.io             import netcdf
import matplotlib.pyplot  as plt
import numpy              as np
dir = '/archive/Junyi.Chai/QG_exp/Aug13/Sc=2/0.1_Fx400_resox2/HISTORY/'
file_set = ['history003.nc','history004.nc']
save_dir = '/home/j1c/analysis/2014/Aug13/movie/Sc=2_drag=0.1/'

NUM_SLICE = 2000
frame     = 0
z_level   = 0 # upper layer = 0; lower layer = 1
C_LIMIT   = 2000
for file_n in range(0,len(file_set)):
    f         = netcdf.netcdf_file(dir+file_set[file_n])
    vor_p     = f.variables['PV_anomaly']
    for t in range(0,NUM_SLICE):
        vor_2d = vor_p[t,z_level,:,:]
        fig = plt.pcolor(vor_2d,cmap=plt.cm.coolwarm)
        plt.clim(-C_LIMIT,C_LIMIT)
        plt.xlim([0, 511])
        plt.ylim([0, 511])
        plt.Axes.set_aspect(plt.axes(),'equal')
        save_name = save_dir + str(frame).zfill(5) + '.png'
        plt.savefig(save_name)
        fig.remove()
        frame = frame + 1
