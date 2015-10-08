#-------------------------------------------------------------------------------
# plot the eddy and total energy spectrum at 3 different times
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/home/j1c/analysis/2015/qg_model/Jan17_drag_1e-4/Ek/'

data1 = np.load(data_dir + 'Ek_seg1.npz')
data2 = np.load(data_dir + 'Ek_seg90.npz')
data3 = np.load(data_dir + 'Ek_seg180.npz')

k = data1['k']

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(k, data1['Ek'], 'b-', label='KE, t1')
ax.loglog(k, data1['EKEk'], 'b--', label='EKE, t1')
ax.loglog(k, data2['Ek'], 'g-', label='KE, t2')
ax.loglog(k, data2['EKEk'], 'g--', label='EKE, t2')
ax.loglog(k, data3['Ek'], 'k-', label='KE, t3')
ax.loglog(k, data3['EKEk'], 'k--', label='EKE, t3')

EKEk  = data3['EKEk']
Ek = data3['Ek']
most_energetic_k = np.where(EKEk == EKEk.max())[0][0]
x   = np.arange(10,101)
p5  = Ek[most_energetic_k]*x**(-5.0)/(x[most_energetic_k]**(-5.0))/5.0
p53 = Ek[most_energetic_k]*x**(-5.0/3.0)/(x[most_energetic_k]**(-5.0/3.0))/100.0
ax.loglog(x,p5,  label='$k^{-5}$')
ax.loglog(x,p53, label='$k^{-5/3}$')
ax.set_xlabel('Wavenumber')
ax.set_ylabel('Energy spectrum')
ax.legend(loc='best')
plt.show()