import nc_tools
import matplotlib.pyplot as plt
import numpy as np
file_Nov4 = ('/archive/Junyi.Chai/QG_exp/Nov4_drag5e-1', 'Nov4_drag5e-1_spectra_seg11.nc')
file_Nov5 = ('/archive/Junyi.Chai/QG_exp/Nov5_Sc1.6_drag5e-1', 'Nov5_Sc1.6_drag5e-1_spectra_seg11.nc')
file_Nov4Reso2x = ('/archive/Junyi.Chai/QG_exp/Nov4_Reso2x_Sc1.6_drag5e-1', 'Nov4_Reso2x_Sc1.6_drag5e-1_spectra_seg6.nc')

var_name = 'xferms'
mode     = 3

var_Nov4 = np.mean(nc_tools.ncread(file_Nov4[0], file_Nov4[1], var_name), 0)
var_Nov5 = np.mean(nc_tools.ncread(file_Nov5[0], file_Nov5[1], var_name), 0)
var_Nov4Reso2x = np.mean(nc_tools.ncread(file_Nov4Reso2x[0], file_Nov4Reso2x[1], var_name), 0)

k = np.arange(1, 512)
k2x = np.arange(1, 1024)
plt.semilogx(k, k*var_Nov4[mode], label='kd=250')
plt.semilogx(k/2, k/2*var_Nov5[mode], label='kd=500, k/2')
plt.semilogx(k2x, k2x*var_Nov4Reso2x[mode], label='kd=250, reso x2')
plt.legend(loc='best')
plt.xlabel('Wavenumber')
#plt.ylabel(r'$U\mathrm{Re}[\psi_{k}^{*}(\partial\nabla^{2}\tau/\partial x)_{k}]$')
plt.ylabel(r'$\mathrm{Re}[\psi_{k}^{*}J_{k}(\tau,\nabla^{2}\tau)]$')
plt.xlim([0, 1023])
plt.show()