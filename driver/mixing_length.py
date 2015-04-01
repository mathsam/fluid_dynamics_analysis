"""
Calculate square_v, heat flux, tracer flux for evaluating mixing length

square_v  = <v, v>, <> denote domain averaged integral
heat_flux = <v, \tau>, \tau is baroclinic stream function 
                       0.5*(phi_upper-phi_lower)
tracer_flux = <v, tracer_y>
"""

import numpy as np
import os
import qg_transform
from nc_tools import NetCDFChain
import matlab_style as ml

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3', '5e-4', '1e-4']
arch_dir  = '/archive/Junyi.Chai/QG_exp/'
save_dir  = '/home/j1c/analysis/2015/qg_model/scales'
last_n_files = 2

square_v    = np.zeros((len(cri_list), len(drag_list)))
square_u    = np.zeros((len(cri_list), len(drag_list)))
square_tau  = np.zeros((len(cri_list), len(drag_list)))
heat_flux   = np.zeros((len(cri_list), len(drag_list)))
tracer_flux = np.zeros((len(cri_list), len(drag_list), 2))


for i, exp_i in enumerate(exp_list):
    for j, drag_i in enumerate(drag_list):
        filedir = arch_dir + exp_i + drag_i
        filename = r'%sseg[0-9]+[0-9]+' %(exp_i + drag_i + '_')
        print filename
        
        psik = NetCDFChain(filedir, filename, 'psi', last_n_files=last_n_files)[:]
        psik = qg_transform.real2complex(psik)
        bt_uk, vk   = qg_transform.get_velocities(psik)
        bt_uk = np.mean(bt_uk, -1)
        bt_vk = np.mean(vk, -1) #barotropic meridional velocity
        square_v[i,j] = np.mean(
            qg_transform.prod_domain_ave_int(bt_vk, bt_vk), 0)

        square_u[i,j] = np.mean(
            qg_transform.prod_domain_ave_int(bt_uk, bt_uk), 0)
        del bt_uk
        
        tau = 0.5*(psik[...,0] - psik[...,1])
        heat_flux[i,j] = np.mean(
            qg_transform.prod_domain_ave_int(bt_vk, tau), 0)
        square_tau[i,j] = np.mean(
            qg_transform.prod_domain_ave_int(tau, tau), 0)
        del tau
        del psik
        del bt_vk
            
        tracer_yk = NetCDFChain(filedir, filename, 'tracer_y', last_n_files=last_n_files)[:]
        tracer_yk = qg_transform.real2complex(tracer_yk)
        tracer_flux[i,j,:] = np.mean(
            qg_transform.prod_domain_ave_int(vk, tracer_yk), 0)

        del tracer_yk, vk

ml.save(os.path.join(save_dir, 'energy_scale.npz'), square_v=square_v,
                                            square_u=square_u,
                                            heat_flux=heat_flux,
                                            tracer_flux=tracer_flux,
                                            square_tau=square_tau) 
