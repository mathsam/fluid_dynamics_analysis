import numpy as np

def expand_field(field2d, d=None):
    """expand the boundary of a doubly periodic field2d
    
    Args:
        field2d: numpy array with shape (n,m)
        d: border size
    Return:
        outfield2d: with shape (n+2*d, m+2*d)
    """
    n, m = field2d.shape
    if d is None:
        d = np.minimum(n, m)
    elif d > np.minimum(n, m):
        raise NotImplementedError("Currently not support to expand to this large")
    outfield2d = np.zeros((n+2*d, m+2*d))
    outfield2d[:d,:d] = field2d[-d:,-d:]
    outfield2d[:d,d:d+m] = field2d[-d:,:]
    outfield2d[:d,d+m:] = field2d[-d:,:d]
    outfield2d[d:d+n,:d] = field2d[:,-d:]
    outfield2d[d:d+n,d:d+m] = field2d
    outfield2d[d:d+n,-d:] = field2d[:,:d]
    outfield2d[d+n:,:d] = field2d[:d,-d:]
    outfield2d[d+n:,d:d+m] = field2d[:d,:]
    outfield2d[d+n:,d+m:] = field2d[:d,:d]
    return outfield2d
    

def extract_local_extreme(field2d, d=10, threashold=0.5):
    """Retain only the local extremas in field2d in the hope to find vortices

    Args:
        field2d: numpy array, assume to be doubly-periodic
        d: define how big the region is to look for local extrema
        threashold: if abs(f)>threashold*abs(f in neighbour region), f is 
                    considered to be local extrema
    Return:
        filtered_field: retain only local extremas; otherwise 0
    """
    n, m = field2d.shape
    filtered_field = np.zeros_like(field2d)
    ext_field2d = expand_field(field2d, d)
    for i in range(0, n):
        for j in range(0, m):
            if field2d[i,j] > threashold*(np.abs(
                ext_field2d[i:i+2*d+1,j:j+2*d+1])).max():
                filtered_field[i,j] = field2d[i,j]
    return filtered_field

def extract_zonal_extreme(field2d, d=10, pctl=95, ratio=2.0):
    """Retain only the extremas compared to the histogram of zonal field
    in the hope to find vortices

    Args:
        field2d: numpy array, assume to be doubly-periodic
        d: define how wide the region is to calculate histogram/percentile
        pctl: field f is retained if 
            abs(f) > ratio*percentile(abs(zonal field), pctl)
        ratio: how far off should be extrema be
    Return:
        filtered_field: retain only extremas; otherwise 0
    """
    n, m = field2d.shape
    ext_field2d = expand_field(field2d, d)
    cutoffs = np.zeros((n,1))
    for i in range(0, n):
        cutoffs[i] = ratio*np.percentile(
            np.abs(ext_field2d[i:i+2*d+1,:].flatten()), pctl)
    retain_index = np.abs(field2d) > cutoffs
    filtered_field = np.zeros_like(field2d)
    filtered_field[:,:] = np.nan
    filtered_field[retain_index] = field2d[retain_index]
    return filtered_field

def cutoff_filter(field2d, abs_threashold=None):
    if abs_threashold is None:
        abs_threashold = field2d.std()
    filtered_field = np.zeros_like(field2d)
    extrema_index = np.abs(field2d)>abs_threashold
    filtered_field[extrema_index] = field2d[extrema_index]
    return filtered_field