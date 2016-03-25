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
ax1 = fig.add_subplot(111)
THREADSHOLD = 0.5
for filename_prefix, label in zip(filename_set, label_set):
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'tracer_spec')
    perct = np.cumsum(-tracer_flux_k)/np.sum(-tracer_flux_k)
    threadshold = np.arange(1, 512)[perct>THREADSHOLD][0]
    ax1.semilogx(k, perct, label=label + ' kc=' + str(threadshold))
    

ax1.plot([1, 511], [THREADSHOLD, THREADSHOLD], '--', label='pct=%.1f' %THREADSHOLD)
ax1.set_ylabel(r"percentage of tracer flux")
ax1.set_xlabel('k')
ax1.legend(loc='best')

plt.show()