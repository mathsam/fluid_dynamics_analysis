import nc_tools
import numpy as np
import matplotlib.pyplot as plt

file_dir = '/archive/Junyi.Chai/sigmaCoord/2014/Apr23_T127_kfe4_corrected'
files = 'Apr23_T127_kfe4_corrected_seg[4-5][0-9]'
ucomp = nc_tools.NetCDFChain(file_dir, files, 'ucomp')
vcomp = nc_tools.NetCDFChain(file_dir, files, 'vcomp')
lat = nc_tools.ncread(file_dir, 'Apr23_T127_kfe4_corrected_seg30.nc', 'lat')

total_times = ucomp.total_time_steps
KE = np.zeros(total_times)
EKE = np.zeros(total_times)
cos_weight = np.cos(np.pi/180*lat)

for t in range(0, total_times):
    print t
    u = ucomp[t,:,:,:]
    v = ucomp[t,:,:,:]
    U2_3d = u**2 + v**2
    U2_1d = np.mean(np.mean(U2_3d, -1), 0)
    KE[t] = np.mean(U2_1d*cos_weight)
    
    umean = np.mean(u, -1)[:,:,np.newaxis]
    vmean = np.mean(v, -1)[:,:,np.newaxis]
    Ueddy2_3d = (u-umean)**2 + (v-vmean)**2
    Ueddy2_1d = np.mean(np.mean(Ueddy2_3d, -1), 0)
    EKE[t] = np.mean(Ueddy2_1d*cos_weight)