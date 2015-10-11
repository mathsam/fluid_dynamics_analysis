import numpy as np
import scipy.ndimage.measurements as measure

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
    filtered_field[:,:] = 0
    filtered_field[retain_index] = field2d[retain_index]
    return filtered_field

def extract_ovals(field2d, min_pixels=9, max_ecc=0.8):
    """retained oval objects in field that has been filtered, e.g., returned 
    from extract_zonal_extreme, which contains disconnected objects.
    
    Args:
        field2d: numpy array returned from filters, e.g. extract_zonal_extreme
        min_pixels: minimum number of pixels a vortex must have
        max_ecc: largest eccentricity a vortex can have
    Return:
        vors: field that contains only vortices, hopefully
    """
    [la, nums] = measure.label(field2d)
    vors = np.zeros_like(field2d)
    for k in range(1,nums):
        x,y = np.where(la==k)
        intensity = field2d[x,y]
        if len(x) >= min_pixels and _eccentricity(x,y,intensity)<=max_ecc:
            vors[x,y] = field2d[x,y]
    return vors

def center_one_vortex(box_size, x, y, with_bkgroud=True, vor=None, *fields):
    """Put a vortex at the center of a box with side length `box_size`
    
    Args:
        box_size: int for size length of the box. Odd number is better
        x: 1d numpy array for x locations
        y: 1d numpy array for y locations
        vor: 2d numpy field of vorticity. vor[x,y] gives vorticity for current
             vortex
             if not present, assume vorticity to be 1 for all (x,y)
        *fields: variable number of other 2d fields, e.g., u, v, psi...
        with_bkgroud: True/False. whether to include points that does not exists
                      in (x,y) pairs.
    Returns:
        if *fields not present:
          cen_vor: 2d numpy array with shape (box_size, box_size), centered 
                   vorticity field
        if *field do present:
          (cen_vor, cen_fld1, cen_fld2, ...)
    """
    if vor is not None:
        momts = _moments2e(x, y, vor[x,y])
        widthx = vor.shape[0]+1
        widthy = vor.shape[1]+1
    else:
        momts = _moments2e(x, y)
        widthx = np.inf
        widthy = np.inf
    vor_ctr_x = int(np.round(momts['mean_x']))
    vor_ctr_y = int(np.round(momts['mean_y']))
    cen_vor = np.zeros((box_size,box_size))
    box_cen = box_size/2
    x_left  = np.maximum(vor_ctr_x-box_cen, 0)
    x_right = np.minimum(vor_ctr_x-box_cen+box_size+1, widthx)
    y_left  = np.maximum(vor_ctr_y-box_cen, 0)
    y_right = np.minimum(vor_ctr_y-box_cen+box_size+1, widthy)
    if vor is not None:
        if with_bkgroud:
            cen_vor = vor[x_left:x_right,
                          y_left:y_right]
        else:
            cen_vor[x-vor_ctr_x+box_cen,y-vor_ctr_y+box_cen] = vor[x,y]
    else:
        cen_vor[x-vor_ctr_x+box_cen,y-vor_ctr_y+box_cen] = 1
    if len(fields) == 0:
        return cen_vor
    opt_fields_out = []
    for i in range(len(fields)):
        curr_field = np.zeros((box_size,box_size))
        if with_bkgroud:
            curr_field = fields[i][x_left:x_right,
                                   y_left:y_right]
        else:
            curr_field[x-vor_ctr_x+box_cen,y-vor_ctr_y+box_cen] = fields[i][x,y]
        opt_fields_out.append(curr_field)
    return (cen_vor,) + tuple(opt_fields_out)
        

def cutoff_filter(field2d, abs_threashold=None):
    if abs_threashold is None:
        abs_threashold = field2d.std()
    filtered_field = np.zeros_like(field2d)
    extrema_index = np.abs(field2d)>abs_threashold
    filtered_field[extrema_index] = field2d[extrema_index]
    return filtered_field

def _eccentricity(x,y, intensity=None):
    """Calculate eccentricity of an image given intensities at (x,y)
    See https://en.wikipedia.org/wiki/Image_moment
    
    Args:
        x: 1d numpy array for x locations
        y: 1d numpy array for y locations
        intensity: intensity of image at location (x,y)
                   if not present, assume intensity to be 1 for all locations
    Return:
        moments: eccentricity of input image. 0 for circle
    """
    moments = _moments2e(x.astype(float),y.astype(float),intensity)
    mup20 = moments['mu20']/moments['mu00']
    mup02 = moments['mu02']/moments['mu00']
    mup11 = moments['mu11']/moments['mu00']
    lambda1 = 0.5*(mup20+mup02 + np.sqrt(4*mup11**2+(mup20-mup02)**2) )
    lambda2 = 0.5*(mup20+mup02 - np.sqrt(4*mup11**2+(mup20-mup02)**2) )
    return np.sqrt(1-lambda2/lambda1)

def _moments2e(x,y, intensity=None):
    """
    This function calculates the raw, centered and normalized moments
    for any image passed as a numpy array.
    
    Further reading:
    https://en.wikipedia.org/wiki/Image_moment
    https://en.wikipedia.org/wiki/Central_moment
    https://en.wikipedia.org/wiki/Moment_(mathematics)
    https://en.wikipedia.org/wiki/Standardized_moment
    http://opencv.willowgarage.com/documentation/cpp/structural_analysis_and_shape_descriptors.html#cv-moments
    
    Args:
        x: 1d numpy array for x locations
        y: 1d numpy array for y locations
        intensity: intensity of image at location (x,y)
                   if not present, assume intensity to be 1 for all locations
    Return:
        moments: dict for up to 3rd moments
    """
    if intensity is None:
        intensity = np.ones_like(x)
    else:
        intensity = np.abs(intensity) #to ensure 2rd moments to be positive
    
    moments = {}
    moments['mean_x'] = np.sum(x*intensity)/np.sum(intensity)
    moments['mean_y'] = np.sum(y*intensity)/np.sum(intensity)
            
    # raw or spatial moments
    moments['m00'] = np.sum(intensity)
    moments['m01'] = np.sum(x*intensity)
    moments['m10'] = np.sum(y*intensity)
    moments['m11'] = np.sum(y*x*intensity)
    moments['m02'] = np.sum(x**2*intensity)
    moments['m20'] = np.sum(y**2*intensity)
    moments['m12'] = np.sum(x*y**2*intensity)
    moments['m21'] = np.sum(x**2*y*intensity)
    moments['m03'] = np.sum(x**3*intensity)
    moments['m30'] = np.sum(y**3*intensity)
    
    # central moments
    moments['mu00'] = moments['m00']
    moments['mu01'] = 0.0
    moments['mu10'] = 0.0
    moments['mu11'] = np.sum((x-moments['mean_x'])*(y-moments['mean_y'])*intensity)
    moments['mu02'] = np.sum((y-moments['mean_y'])**2*intensity) # variance
    moments['mu20'] = np.sum((x-moments['mean_x'])**2*intensity) # variance
    moments['mu12'] = np.sum((x-moments['mean_x'])*(y-moments['mean_y'])**2*intensity)
    moments['mu21'] = np.sum((x-moments['mean_x'])**2*(y-moments['mean_y'])*intensity) 
    moments['mu03'] = np.sum((y-moments['mean_y'])**3*intensity) 
    moments['mu30'] = np.sum((x-moments['mean_x'])**3*intensity) 
    
    
    # opencv versions
    #moments['mu02'] = np.sum(intensity*(x-m01/m00)**2)
    #moments['mu02'] = np.sum(intensity*(x-y)**2)
    
    # wiki variations
    #moments['mu02'] = m20 - mean_y*m10 
    #moments['mu20'] = m02 - mean_x*m01
        
    # central standardized or normalized or scale invariant moments
    moments['nu11'] = moments['mu11'] / np.sum(intensity)**(2/2+1)
    moments['nu12'] = moments['mu12'] / np.sum(intensity)**(3/2+1)
    moments['nu21'] = moments['mu21'] / np.sum(intensity)**(3/2+1)
    moments['nu20'] = moments['mu20'] / np.sum(intensity)**(2/2+1)
    moments['nu03'] = moments['mu03'] / np.sum(intensity)**(3/2+1) # skewness
    moments['nu30'] = moments['mu30'] / np.sum(intensity)**(3/2+1) # skewness
    return moments