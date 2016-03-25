import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab

filename_set = ['Dec12_kf16_drag256e-4',
    'Dec12_kf32_drag256e-4',
    'Dec12_kf64_drag256e-4',
    'Dec12_kf128_drag256e-4']

# filename_set = ['Dec12_kf16_drag8e-4',
#     'Dec12_kf32_drag8e-4',
#     'Dec12_kf64_drag8e-4']
label_set = ['kf=16','kf=32','kf=64','kf=128']
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
for filename_prefix, label in zip(filename_set, label_set):
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'ty_spec')
    mlab.load(save_dir + 'tracer_spec')
    ax1.loglog(k, -tracer_flux_k, label=label)
    ax2.loglog(k, -tracer_flux_k/np.sqrt(v2_k*eddy_t2[:,np.newaxis]), label=label)

ax1.set_ylabel(r"$v'S'$ spec")
ax1.legend(loc='best')
ax2.set_ylabel(r"corr($v'$,$S'$)")
ax2.legend(loc='best')
plt.show()