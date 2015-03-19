import numpy as np

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