import qg_transform
import nc_tools
import numpy as np
import diag

# Sc = 1.6, F = 1583.14, beta = 3957.86
# Sc = 2.0, F = 1583.14, beta = 3166.29
# Sc = 2.5, F = 1583.14, beta = 2533.03
F    = 1583.14
beta = 3957.86

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=1)[:]
psi = qg_transform.real2complex(psi)
v = qg_transform.get_velocities(psi)[1]
v = np.mean(v, -1)[..., np.newaxis]

##
v_large = qg_transform.filter_zonal(v, 25, None, True)
v_small = qg_transform.filter_zonal(v, None, 24, True)
v_largeg = qg_transform.spec2grid(v_large)
v_smallg = qg_transform.spec2grid(v_small)
mean_v2l = np.mean(np.mean(v_largeg**2, 0), -2)
mean_v2s = np.mean(np.mean(v_smallg**2, 0), -2)

##
climitology = np.load('/home/j1c/analysis/2015/qg_model/Jan17_drag_1e-4/zonal_mean_field.npz')

##
import matlab_style as mlab
import matplotlib.pyplot as plt
indx_s=420
indx_e=700

bt_U = np.mean(climitology['mean_u'], -1)
fig = plt.figure()
ax1 = fig.add_subplot(111)
y = np.linspace(-np.pi, np.pi, 1024)
ax1_twin = mlab.plotyy(y[indx_s:indx_e], bt_U[indx_s:indx_e], y[indx_s:indx_e], mean_v2l[indx_s:indx_e], None, r'$U_{bt}$', r"$v_{bt}'^2$, kx>=25", ax1)
ax1_twin.plot(y[indx_s:indx_e], mean_v2s[indx_s:indx_e], label=r"$v_{bt}'^2$, kx<=24")
ax1.legend(loc='upper left')
ax1_twin.legend(loc='upper right', fancybox=False)
ax1_twin.yaxis.grid(True)
plt.show()