import nc_tools
import matplotlib.pyplot as plt
import numpy as np


archive_dir = '/archive/Junyi.Chai/shallow_water/'
exp_prefix  = 'Oct17_amp'
exp_list = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9',
            '1','2','3','4','5','6','7','8','9','10',
            '11','12','13','14','15','16','17','18','19','20',
            '21','22','23','24','25','26','27','28','29','30']

tr_flux_at45 = np.zeros(len(exp_list))
tr_flux_NH = np.zeros(len(exp_list))
Vstd_NH = np.zeros(len(exp_list))

for i in range(len(exp_list)):
    exp_name = exp_prefix + str(exp_list[i])
    file_dir = archive_dir + exp_name + '/'
    nc_file  =  exp_name + '_seg1.nc'
    vcomp = nc_tools.ncread(file_dir, nc_file, 'vcomp')
    vcomp_mean = np.mean(vcomp, -1)[:,:,np.newaxis]
    vcomp -= vcomp_mean
    tr    = nc_tools.ncread(file_dir, nc_file, 'tr')
    v_tr  = np.mean(np.mean(tr*vcomp, -1), 0)
    EKE1d = np.mean(vcomp[0]**2, -1)
    EKE1d_weighted = EKE1d*cos_lat
    Vstd_NH[i] = np.sqrt(np.mean(EKE1d_weighted[64:]))
    #plt.plot(lat, v_tr)
    tr_flux_at45[i] = v_tr[96]
    tr_flux_NH[i] = np.mean(v_tr[64:])
    
##
file_dir = '/archive/Junyi.Chai/shallow_water/Oct17_amp1'
nc_file = 'Oct17_amp1_seg1.nc'
vcomp = nc_tools.ncread(file_dir, nc_file, 'vcomp')
vcomp_mean = np.mean(vcomp, -1)[:,:,np.newaxis]
vcomp -= vcomp_mean
tr    = nc_tools.ncread(file_dir, nc_file, 'tr')
tr_mean = np.mean(tr, -1)[:,:,np.newaxis]
tr -= tr_mean
v_tr  = np.mean(np.mean(tr*vcomp, -1), 0)
#plt.plot(lat[64:], v_tr[64:])
plt.plot(lat[64:], v_tr[64:]/v_tr_base[64:])
#plt.show()

    
## amplitude vs. flux
plt.loglog(x, -tr_flux_NH,label='wave amplitide vs. eddy tracer flux')
plt.loglog(x, 0.0019051874548182781*x, '--', label='y ~ x')
plt.loglog(x, 0.0003012695802963167*x**2, '--', label='y ~ x^2')

plt.xlabel('Wave amplitude')
plt.ylabel('Eddy tracer flux')
plt.legend(loc='best')
plt.show()
## y structure of disturbances
vor = nc_tools.ncread('/archive/Junyi.Chai/shallow_water/Oct17_amp1', 'Oct17_amp1_seg1.nc', 'tr')
lat = nc_tools.ncread('/archive/Junyi.Chai/shallow_water/Oct17_amp1', 'Oct17_amp1_seg1.nc', 'lat')
zonalmean_vor = np.mean(vor, -1)[:,:,np.newaxis]
eddy_vor = vor-zonalmean_vor
plt.contourf(lat[64:], range(15), eddy_vor[:15,64:,1])
plt.draw()

##
file_dir = '/archive/Junyi.Chai/shallow_water/Oct17_amp10'
nc_file = 'Oct17_amp10_seg1.nc'
vcomp = nc_tools.ncread(file_dir, nc_file, 'vcomp')
vcomp_mean = np.mean(vcomp, -1)[:,:,np.newaxis]
vcomp -= vcomp_mean
tr    = nc_tools.ncread(file_dir, nc_file, 'tr')
v_tr  = np.mean(np.mean(tr*vcomp, -1), 0)