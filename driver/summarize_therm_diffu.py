## summarize energy generation rate and evaulate heat diffusivity
from pandas import Series, DataFrame
import pandas as pd
from math import pi
from nc_tools import NetCDFChain

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3']
arch_dir  = '/archive/Junyi.Chai/QG_exp/'
average_period = 200
U_BAR    = 1.
ROSSBY_R = 2*pi/250


gen_rate = DataFrame(columns=cri_list, index=map(float,drag_list))
gen_rate.index.name   = 'bottom_drag'
gen_rate.columns.name = 'criticality'

for i, exp_i in enumerate(exp_list):
    for drag_i in drag_list:
        filedir = arch_dir + exp_i + drag_i
        filename = r'%senergy_seg[0-9]+' %(exp_i + drag_i + '_')
        gen_bci_rate     = NetCDFChain(filedir, filename, 'gen_bci_rate')[-200:]
        ave_gen_bci_rate = gen_bci_rate.mean()
        gen_rate[cri_list[i]][float(drag_i)] = ave_gen_bci_rate

therm_diffusivity = gen_rate/(2*pi)**2*ROSSBY_R/U_BAR**3

## save results
import os
save_dir = '/home/j1c/analysis/2015/qg_model/Jan17to18_summary/'
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

therm_diffusivity.to_csv(save_dir + 'therm_diffusivity.csv')
therm_diffusivity.save(save_dir + 'therm_diffusivity.pd')
gen_rate.to_csv(save_dir + 'gen_rate.csv')
gen_rate.save(save_dir + 'gen_rate.pd')
## plot results
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
therm_diffusivity.plot(ax=ax, style='--o', loglog=True)
ax.set_ylabel('thermal diffusivity')
plt.show()
fig.savefig(save_dir + 'diffusivity.png')