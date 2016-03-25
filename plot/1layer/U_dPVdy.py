import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab

filename_set = ['Dec12_kf16_drag256e-4',
    'Dec12_kf32_drag256e-4',
    'Dec12_kf64_drag256e-4']
# label_set = ['kf=16','kf=32','kf=64','kf=128']
# filename_set = ['Dec12_kf16_drag8e-4',
#     'Dec12_kf32_drag8e-4',
#     'Dec12_kf64_drag8e-4']
beta = 50.2654824574
num_ks = len(filename_set)
y = np.linspace(-np.pi, np.pi, 1024)
fig = plt.figure()

for i, filename_prefix in enumerate(filename_set):
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'noboru_budget')
    ax1 = fig.add_subplot(num_ks,1, i+1)
    ax1.plot(y, dfdy(lang_btpv), label='noboru dPVdy')
    ax1.set_ylim([0, 150])
    ax2 = mlab.plotyy(y, -dfdy(dfdy(ug1d))+beta, y, ug1d, ylabel1='dPVdy', ylabel2='U', ax1=ax1)
    ax1.legend(loc='upper left', fontsize=9)
    ax2.legend(loc='best', fontsize=9)
    ax2.set_ylim([-2, 2])

plt.show()
