## summarize energy generation rate and evaulate heat diffusivity
from pandas import Series, DataFrame
import pandas as pd
from math import pi
from nc_tools import NetCDFChain

exp_list  = ['Nov4_drag', 'Nov5_Sc1.6_drag']
cri_list  = ['1.6', '1.6']
drag_list = ['5e-1', '5e-2', '5e-3', '5e-4', '5e-5']
arch_dir  = '/archive/Junyi.Chai/QG_exp/'
average_period = 50
U_BAR    = 1.
ROSSBY_R = 2*pi/500.

num_rows = len(drag_list)
num_colmns = len(cri_list)

energetics = {'ke': np.zeros((num_rows, num_colmns)),
    'filter_rate': np.zeros((num_rows, num_colmns)),
    'bottom_drag_rate': np.zeros((num_rows, num_colmns)),
    'ape': np.zeros((num_rows, num_colmns)),
    'gen_bci_rate': np.zeros((num_rows, num_colmns)),
    'eddy_time': np.zeros((num_rows, num_colmns))}


for exp_i, exp_name in enumerate(exp_list):
    for drag_i, drag_name in enumerate(drag_list):
        filedir = arch_dir + exp_name + drag_name
        filename = r'%senergy_seg[0-9]+' %(exp_name + drag_name + '_')
        for var_name in energetics.keys():
            var = NetCDFChain(filedir, filename, var_name)[-average_period:]
            ave_var = np.mean(var)
            energetics[var_name][drag_i, exp_i] = ave_var

#therm_diffusivity = gen_rate/(2*pi)**2*ROSSBY_R/U_BAR**3

## save results
import os
save_dir = '/home/j1c/analysis/2015/qg_model/Nov4_5/'
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

for var_name in energetics.keys():
    np.save(save_dir + var_name, energetics[var_name])
## plot result
plt.loglog(drag, energetics['gen_bci_rate'][:,1], '-ro', label='energy generation rate')
plt.loglog(drag, 100*drag, '--', label='~ x')
plt.loglog(drag, 1000*np.sqrt(drag), '--', label='~ x^{1/2}')
plt.legend(loc='best')
plt.xlabel('drag')
plt.ylabel(r'$\epsilon$')
plt.show()