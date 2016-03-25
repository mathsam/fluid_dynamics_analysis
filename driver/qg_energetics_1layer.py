## load data
import sys
sys.path.append('/home/j1c/py/lib');
#import matlab_style as ml
#ml.clear()
from nc_tools import NetCDFChain
import numpy as np

filename_prefix = 'Mar14_kf128_drag4096e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/2016/%s' %filename_prefix
filename = r'%s_energy_seg[0-9]+' %filename_prefix


t   = NetCDFChain(filedir, filename, 'time')[:]
ke  = NetCDFChain(filedir, filename, 'ke')[:]
gen_rmf_rate     = NetCDFChain(filedir, filename, 'gen_rmf_rate')[:]
bottom_drag_rate = NetCDFChain(filedir, filename, 'bottom_drag_rate')[:]
filter_rate      = NetCDFChain(filedir, filename, 'filter_rate')[:]

mean_gen = np.mean(gen_rmf_rate[-100:])
mean_fric = np.mean(bottom_drag_rate[-100:])
mean_diss= np.mean(bottom_drag_rate[-100:]+filter_rate[-100:])
print "gen = %f, diss=%f, imblance=%f" %(mean_gen, mean_diss, (mean_gen+mean_diss)/mean_gen)
print "frictional dissipation rate=%.3e" %-mean_fric
# plotting
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2,1)
axes[0].plot(t, ke,  label='ke')
axes[0].set_xlabel('time')
axes[0].set_ylabel('energy')
axes[0].legend(loc='best')
    
axes[1].plot(t, gen_rmf_rate, label='gen rate')
axes[1].plot(t, -bottom_drag_rate, label='- drag rate')
axes[1].plot(t, -filter_rate, label='- filter rate')
axes[1].set_xlabel('time')
axes[1].set_ylabel('energy')
axes[1].legend(loc='best')
plt.show()

## tracer_y
filename_prefix = 'Dec12_kf128_drag64e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
filename = r'%s_energy_seg([0-9]+)' %filename_prefix

t         = NetCDFChain(filedir, filename, 'time')[:]
tvary     = NetCDFChain(filedir, filename, 'tvary')[:]
gen_ty    = NetCDFChain(filedir, filename, 'gen_tg_ty')[:]
filter_ty = NetCDFChain(filedir, filename, 'filter_rate_ty')[:]

gen_ty_mean    = np.mean(gen_ty[-100:])
filter_ty_mean = np.mean(filter_ty[-100:])
print "gen_ty=%.2e, filter_ty=%.2e, imbalance=%f)" %(gen_ty_mean, filter_ty_mean,
                                    (gen_ty_mean+filter_ty_mean)/gen_ty_mean)
fig_t, axes_t = plt.subplots(2,1)
axes_t[0].plot(t, tvary,  label='tracer variance')
axes_t[0].set_xlabel('time')
axes_t[0].set_ylabel('tracer variance')
axes_t[0].legend(loc='best')

axes_t[1].plot(t, gen_ty, label='gen tracer var')
axes_t[1].plot(t, -filter_ty, label='- filter')
axes_t[1].set_xlabel('time')
axes_t[1].set_ylabel('dvar/dt')
axes_t[1].legend(loc='best')
plt.show()                        

## save figure
import os
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
    
save_name = 'energetics.png'
fig.savefig(save_dir + save_name)


