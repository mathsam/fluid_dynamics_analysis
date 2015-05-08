import numpy as np
import qg_transform

def centroid(k, weight):
    """
    calculate the centroid of an array k given its weight, defined as
        ave_k = mean(k*weight)/mean(weight)
    
    Args:
        k: 1d numpy array
        weight: same dimension as k, non-negative numbers
        
    Returns:
        ave_k: centroid of k
    """
    if k.shape != weight.shape:
        raise TypeError('k and weight must have the same type')
    ave_k = np.mean(k*weight)/np.mean(weight)
    return ave_k
    
def inverse_centroid(k, weight, order=1.0):
    """
    giving more weights to smaller k, defined as
        ave_k^(-order) = mean(k^(-order) * weight) / mean(weight)
        
    Args:
        k: 1d numpy array
        weight: same dimension as k, non-negative numbers
        
    Returns:
        ave_k: inverse centroid of k
    """
    if k.shape != weight.shape:
        raise TypeError('k and weight must have the same type')
    order = float(order)
    i_ave_k = np.mean(k**(-order) * weight)/np.mean(weight)
    ave_k   = i_ave_k**(-1./order)
    return ave_k
    
# zonal mean zonal wind                      
def zonal_mean_zonal_wind(psi):
    timeave_psi = np.mean(psi, 0)
    uk, vk = qg_transform.get_velocities(timeave_psi)
    ug  = qg_transform.spec2grid(uk)    
    mean_ug = np.mean(ug, -2)
    return mean_ug

# zonal mean PV
def zonal_mean_PV(psi, F, beta):
    timeave_psi = np.mean(psi, 0)
    pvg = qg_transform.spec2grid(qg_transform.get_PV(timeave_psi, F))
    pvg += qg_transform.get_betay(pvg, beta)    
    mean_pvg = np.mean(pvg, -2)
    return mean_pvg

def zonal_mean_dPVdy(psi, F, beta):
    timeave_psi = np.mean(psi, 0)
    pvk = qg_transform.get_PV(timeave_psi, F)
    dpvkdy = qg_transform.partial_y(pvk)
    dpvgdy = qg_transform.spec2grid(dpvkdy)
    mean_dpvdy = np.mean(dpvgdy, -2) + beta
    return mean_dpvdy
    
def sponge_mask(num_lats, sponge_width=0.05, sponge_rate=1.):
    """return the sponge_mask for calculating energy dissipation by the sponges
    Args:
        sponge_width: ratio of the width of the sponge to the domain size
        num_lats: number of grid points
        sponge_rate: default is 1
    Returns:
        mask: numpy array with shape (num_lats, )
    """
    lat = np.linspace(0, 1, num_lats)
    lower_sponge = sponge_rate*(1. - lat/sponge_width)
    upper_sponge = sponge_rate*(lat - 1. + sponge_width)/sponge_width
    
    mask = np.maximum(lower_sponge, 
        np.maximum(upper_sponge, np.zeros(num_lats)))
    return mask