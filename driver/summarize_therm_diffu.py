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


gen_rate = DataFrame(columns=cri_list, index=drag_list)

for i, exp_i in enumerate(exp_list):
    for drag_i in drag_list:
        filedir = arch_dir + exp_i + drag_i
        filename = r'%senergy_seg[0-9]+' %(exp_i + drag_i + '_')
        gen_bci_rate     = NetCDFChain(filedir, filename, 'gen_bci_rate')[-200:]
        ave_gen_bci_rate = gen_bci_rate.mean()
        gen_rate[cri_list[i]][drag_i] = ave_gen_bci_rate

therm_diffusivity = gen_rate/(2*pi)**2*ROSSBY_R/U_BAR**3