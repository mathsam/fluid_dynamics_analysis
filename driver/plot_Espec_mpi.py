from scipy.io import netcdf
from numpy    import *
import matplotlib.pyplot as plt

filename = '/archive/Junyi.Chai/QG_exp/Jan17_drag_1e-3/Jan17_drag_1e-3_seg40.nc'

f = netcdf.netcdf_file(filename,'r',mmap=False)
psi_spec_p = f.variables['psi']

T_START = 1  # starting from 1
T_END   = psi_spec_p.shape[0]
z_level = 0

nkx      = psi_spec_p.shape[3]
nky      = psi_spec_p.shape[2]
kmax     = nky-1
indx_kx0 = (nkx-1)/2
field    = zeros((nky,nkx))
Ek       = zeros((kmax,1))
EkEKE    = zeros((kmax,1))
ksqd_    = zeros((nky,nkx))
fieldEKE = zeros((nky,nkx))

for i in range(0,nkx):
    for j in range(0,nky):
        ksqd_[j][i] = (i-kmax)**2 + j**2
        
radius_arr = floor(sqrt(ksqd_))

for i in range(T_START-1,T_END):
    field = field + ksqd_*(psi_spec_p[i,0,:,:,z_level]**2+psi_spec_p[i,1,:,:,z_level]**2)

field    = field/T_END
fieldEKE = copy(field)
fieldEKE[:,indx_kx0] = 0

for i in range(0,kmax):
    Ek[i]     = sum(field[radius_arr == i+1])
    EkEKE[i]  = sum(fieldEKE[radius_arr == i+1])

k_set = arange(1,kmax+1)

## ploting
plt.loglog(k_set,Ek)
plt.loglog(k_set[1:],EkEKE[1:],ls='--')
x   = arange(10,101)
p5  = Ek[9]*x**(-5.0)/(x[0]**(-5.0))/5.0
p53 = Ek[9]*x**(-5.0/3.0)/(x[0]**(-5.0/3.0))*5.0
plt.loglog(x,p5)
plt.loglog(x,p53)
plt.xlim([1,300])
plt.show()