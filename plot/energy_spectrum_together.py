import os
import numpy as np
import diag
import matplotlib.pyplot as plt
import matlab_style as ml
analysis_dir = '/home/j1c/analysis/2015/qg_model'
save_dir = '/home/j1c/analysis/2015/qg_model/heat_flux'

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3', '5e-4', '1e-4']
nondim_drag = ['1.1e-1', '1.1e-2', '1.1e-3', '1.1e-4', '5.6e-5', '1.1e-5']
num_wavenumbers = 511

plt.figure(figsize=(8,8))
for i in range(1,2): # Sc
    for j in [1, 2, 5]: # friction
        exp_name = exp_list[i] + drag_list[j]
        file = os.path.join(analysis_dir, exp_name, 'Ek.npz')
        spec_zip = np.load(file)        
        
        k = spec_zip['k']
        Ek = spec_zip['Ek']
        plt.loglog(k, Ek, label='Sc=' + cri_list[i] 
            + r', $r\lambda /U=$' + nondim_drag[j])
            
#
most_energetic_k = 10
x   = np.arange(10,101)
xs = np.arange(50, 300)
p5  = Ek[most_energetic_k]*x**(-5.0)/(x[most_energetic_k]**(-5.0))/40
p3  = Ek[most_energetic_k]*xs**(-3.0)/(xs[most_energetic_k]**(-3.0))/40
plt.loglog(x,p5,  label='$k^{-5}$')
plt.loglog(xs,p3,  label='$k^{-3}$')
plt.legend(loc='best')
plt.grid()
plt.xlabel('Wavenumber')
plt.ylabel('E(k)')
plt.show()

## EKE spectrum
plt.figure(figsize=(8,8))
for i in range(0,1): # Sc
    for j in [1, 2, 5]: # friction
        exp_name = exp_list[i] + drag_list[j]
        file = os.path.join(analysis_dir, exp_name, 'Ek.npz')
        spec_zip = np.load(file)        
        
        k = spec_zip['k']
        Ek = spec_zip['EKEk']
        plt.loglog(k, Ek, label='Sc=' + cri_list[i] 
            + r', $r\lambda /U=$' + nondim_drag[j])
            
#
most_energetic_k = 10
x   = np.arange(10,101)
p53  = Ek[most_energetic_k]*x**(-5.0/3)/(x[most_energetic_k]**(-5.0/3))
plt.loglog(x,p53,  label='$k^{-5/3}$')
plt.legend(loc='best')
plt.grid()
plt.xlabel('Wavenumber')
plt.ylabel('E(k)')
plt.show()