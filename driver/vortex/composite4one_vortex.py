import vortex_track.image_process as img_proc
import scipy.ndimage.measurements as measure
import operator

box_size = 15
fields_list = [vorg[:,:,1], vorg[:,:,0], 
               ug[:,:,1], ug[:,:,0], 
               vg[:,:,1], vg[:,:,0]]
p_vortex = []
n_vortex = []
num_p_vor = 0
num_n_vor = 0
num_fields = len(fields_list)
for i in range(num_fields):
    p_vortex.append(np.zeros((box_size, box_size)))
    n_vortex.append(np.zeros((box_size, box_size)))

filtered_vor = img_proc.extract_zonal_extreme(vor, 10, 80, 2.0)
vortices = img_proc.extract_ovals(filtered_vor, 9, 0.8)
[la, nums] = measure.label(vortices)
for k in range(0, nums):
    x, y = np.where(la==k)
    vor_in_box = img_proc.center_one_vortex(box_size, x, y, True, *fields_list)
    if np.max(vor_in_box) > 500:
        num_p_vor += 1
        for i in range(num_fields):
            p_vortex[i] += vor_in_box[i]
    elif np.min(vor_in_box) < -500:
        num_n_vor += 1
        for i in range(num_fields):
            n_vortex[i] += vor_in_box[i]
        
p_one_vortex = map(lambda x: x/num_p_vor, p_vortex)
n_one_vortex = map(lambda x: x/num_n_vor, n_vortex)

##
plt.subplot(1,2,1)
plt.imshow(p_one_vortex[0], interpolation='none',cmap='gray')
plt.title('lower')
plt.colorbar()
plt.subplot(1,2,2)
plt.imshow(p_one_vortex[1], interpolation='none',cmap='gray')
plt.title('upper')
plt.colorbar()
plt.show()