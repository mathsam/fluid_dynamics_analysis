## load data
import sys
sys.path.append('/home/j1c/py/lib');
#import matlab_style as ml
#ml.clear()
from nc_tools import NetCDFChain

filename_prefix = 'Nov4_Sc2.5_drag5e-2'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_energy_seg[0-9]+' %filename_prefix


t   = NetCDFChain(filedir, filename, 'time')[:]
ke  = NetCDFChain(filedir, filename, 'ke')[:]
ape = NetCDFChain(filedir, filename, 'ape')[:]
gen_bci_rate     = NetCDFChain(filedir, filename, 'gen_bci_rate')[:]
bottom_drag_rate = NetCDFChain(filedir, filename, 'bottom_drag_rate')[:]
filter_rate      = NetCDFChain(filedir, filename, 'filter_rate')[:]


# plotting
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2,1)
axes[0].plot(t, ke,  label='ke')
axes[0].plot(t, ape, label='ape')
axes[0].set_xlabel('time')
axes[0].set_ylabel('energy')
axes[0].legend(loc='best')
    
axes[1].plot(t, gen_bci_rate, label='gen rate')
axes[1].plot(t, -bottom_drag_rate, label='- drag rate')
axes[1].plot(t, -filter_rate, label='- filter rate')
axes[1].set_xlabel('time')
axes[1].set_ylabel('energy')
axes[1].legend(loc='best')
plt.show()

## save figure
import os
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
    
save_name = 'energetics.png'
fig.savefig(save_dir + save_name)


