import os
import sys
sys.path.append('../lib')
import nc_tools 
import animation_tools
import matplotlib.animation as animation
#from matplotlib import cm

exp_name = 'Mar25_amp5e-2'
z_level  = 0
dpi      = 900
save_dir  = '/home/j1c/analysis/2016/qg_model/%s' %exp_name
save_file = save_dir + '/vor_perb.mp4' 
psi=nc_tools.ncread('/archive/Junyi.Chai/QG_exp/2015/%s' %exp_name, 
                    '%s_seg[0-9]+.nc' %exp_name, 'psi')
frames  = animation_tools.AnimationFrames(psi, z_level,True)
fig,ani = animation_tools.make_animation(frames, frames.total_frames)
## set encoder and dpi
writer = animation.writers['ffmpeg'](fps=25)
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
if not os.path.exists(save_file):
    ani.save(save_file, writer=writer, dpi=dpi)
else:
    raise IOError('file already exists')
