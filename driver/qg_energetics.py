## plotting
def plot_energetics():
    import matplotlib.pyplot as plt
    plt.subplot(2,1,1)
    plt.plot(t, ke,  label='ke')
    plt.plot(t, ape, label='ape')
    plt.xlabel('time')
    plt.ylabel('energy')
    plt.legend(loc='best')
    
    plt.subplot(2,1,2)
    plt.plot(t, gen_bci_rate, label='gen rate')
    plt.plot(t, -bottom_drag_rate, label='- drag rate')
    plt.plot(t, -filter_rate, label='- filter rate')
    plt.xlabel('time')
    plt.ylabel('energy')
    plt.legend(loc='best')
    #plt.show()
    plt.savefig(save_dir + '/' + save_name)
    plt.close()

## load data
import sys
sys.path.append('/home/j1c/py/lib');
#import matlab_style as ml
#ml.clear()
from nc_tools import nc_chains

filedir  = '/archive/Junyi.Chai/QG_exp/Nov28_kfe2_qg'
filename = r'Nov28_kfe2_qg_energy_seg[0-9]+'

save_dir  = '/home/j1c/analysis/2014/QG_model/Nov28_ke2'
save_name = 'energetics_seg100.png'

t   = nc_chains(filedir, filename, 'time')[:]
ke  = nc_chains(filedir, filename, 'ke')[:]
ape = nc_chains(filedir, filename, 'ape')[:]
gen_bci_rate     = nc_chains(filedir, filename, 'gen_bci_rate')[:]
bottom_drag_rate = nc_chains(filedir, filename, 'bottom_drag_rate')[:]
filter_rate      = nc_chains(filedir, filename, 'filter_rate')[:]

plot_energetics()



