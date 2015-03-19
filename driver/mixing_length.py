"""
Calculate std_v, heat flux, tracer flux for evaluating mixing length

std_v     = <v, v>, <> denote domain averaged integral
heat_flux = <v, \tau>, \tau is baroclinic stream function 
                       0.5*(phi_upper-phi_lower)
tracer_flux = <v, tracer_y>
"""

import numpy as np
import qg_transform
from nc_tools import NetCDFChain

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3']
arch_dir  = '/archive/Junyi.Chai/QG_exp/'
AVE_PERIOD = 200
U_BAR    = 1.
ROSSBY_R = 2*np.pi/250


std_v       = np.zeros((len(cri_list), len(drag_list)))
heat_flux   = np.zeros((len(cri_list), len(drag_list)))
tracer_flux = np.zeros((len(cri_list), len(drag_list), 2))

def cal_time_ave_std(prod_int):
    return np.std(np.mean(prod_int, 0))

for i, exp_i in enumerate(exp_list):
    for j, drag_i in enumerate(drag_list):
        filedir = arch_dir + exp_i + drag_i
        filename = r'%sseg[0-9]+' %(exp_i + drag_i + '_')
        print filename
        
        psik = NetCDFChain(filedir, filename, 'psi')[-AVE_PERIOD:]
        psik = qg_transform.real2complex(psik)
        _, vk   = qg_transform.get_velocities(psik)
        bt_vk = np.mean(vk, -1) #barotropic meridional velocity
        std_v[i,j] = cal_time_ave_std(
            qg_transform.prod_domain_ave_int(bt_vk, bt_vk))
        
        tau = 0.5*(psik[...,0] - psik[...,1])
        heat_flux[i,j] = cal_time_ave_std(
            qg_transform.prod_domain_ave_int(bt_vk, tau))
            
        tracer_yk = NetCDFChain(filedir, filename, 'tracer_y')[-AVE_PERIOD:]
        tracer_yk = qg_transform.real2complex(tracer_yk)
        tracer_flux[i,j,0] = cal_time_ave_std(
            qg_transform.prod_domain_ave_int(vk[...,0], tracer_yk[...,0]))
        tracer_flux[i,j,1] = cal_time_ave_std(
            qg_transform.prod_domain_ave_int(vk[...,1], tracer_yk[...,1]))
            
## calculate eddy scale and jet scale
# ------------------------------------------------------------------------------
# eddy scale is defined as (inverse) centroid of EKE
# jet scale is defined as (inverse) centroid of ZKE
# ------------------------------------------------------------------------------
import os
import diag
analysis_dir = '/home/j1c/analysis/2015/qg_model'

cent_EKE = np.zeros((len(cri_list), len(drag_list)))
cent_ZKE = np.zeros((len(cri_list), len(drag_list)))

for i, exp_i in enumerate(exp_list):
    for j, drag_i in enumerate(drag_list):
        exp_name = exp_i + drag_i
        file_Ek = os.path.join(analysis_dir, exp_name, 'Ek.npz')
        Ekzip = np.load(file_Ek)
        k    = Ekzip['k']
        EKEk = Ekzip['EKEk']
        Ek   = Ekzip['Ek']
        ZKEk = Ek - EKEk
        
        cent_EKE[i,j] = diag.centroid(k, EKEk)
        cent_ZKE[i,j] = diag.centroid(k, ZKEk)