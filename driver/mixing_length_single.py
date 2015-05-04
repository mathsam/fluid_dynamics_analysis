import numpy as np
import os
import qg_transform
import nc_tools
import matlab_style as ml

filename_prefix = 'Jan18_c2.5_drag_1e-3'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = '%s_seg[0-9]+[0-9]+' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix

psik = nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=2)[:]
psik = qg_transform.real2complex(psik)
bt_uk, vk   = qg_transform.get_velocities(psik)
bt_vk = np.mean(vk, -1) #barotropic meridional velocity
square_v = np.mean(
    qg_transform.prod_domain_ave_int(bt_vk, bt_vk), -1)


tau = 0.5*(psik[...,0] - psik[...,1])
heat_flux = np.mean(
    qg_transform.prod_domain_ave_int(bt_vk, tau), -1)

#del tau
#del psik
#del bt_vk
    
tracer_yk = nc_tools.ncread(filedir, filename,'tracer_y')[:]
tracer_yk = qg_transform.real2complex(tracer_yk)

tracer_flux = np.mean(
    qg_transform.prod_domain_ave_int(vk, tracer_yk), 0)

#del tracer_yk, vk
##
#save npz
import matlab_style as mlab
import os

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

mlab.save(save_dir + 'mixing_length', square_v=square_v, heat_flux=heat_flux, 
                                      tracer_flux=tracer_flux)