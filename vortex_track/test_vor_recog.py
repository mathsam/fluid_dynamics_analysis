import scipy.ndimage.measurements as measure
f_vorg = extract_zonal_extreme(vorg[:,:,1],10,80)
[la,nums] = measure.label(f_vorg)
for k in range(1,nums):
    x,y = np.where(la==k)
    field2d = vorg[:,:,1]
    intensity = field2d[x,y]
    if len(x) >= 9 and _eccentricity(x,y,intensity) < 0.8:
        plt.subplot(1,2,1)
        plt.scatter(y,x,np.abs(intensity))
        plt.gca().invert_yaxis()
        plt.title("ecc cont=%f, bin=%f" %(_eccentricity(x,y,intensity),
                                                _eccentricity(x,y)))
        plt.subplot(1,2,2)
        x_min = np.maximum(1, np.min(x)-10)
        x_max = np.minimum(field2d.shape[0],np.max(x)+10)
        y_min = np.maximum(1, np.min(y)-10)
        y_max = np.minimum(field2d.shape[1],np.max(y)+10)
        plt.imshow(field2d[x_min:x_max,y_min:y_max], interpolation='none', cmap='gray')
        plt.colorbar()
        plt.scatter(y-y_min,x-x_min)
        plt.savefig('%d.png' %k)
        plt.clf()