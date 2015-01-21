import nc_tools
import animation_tools
z_level = 0
psi=nc_tools.ncread('/archive/Junyi.Chai/QG_exp/Jan17_drag_1e-3','Jan17_drag_1e-3_seg[0-9]+','psi')
frames  = animation_tools.AnimationFrames(psi, z_level)
fig,ani = animation_tools.make_animation(frames, frames.total_frames)
##
import matplotlib.animation as animation
writer = animation.writers['ffmpeg'](fps=10)
ani.save('/home/j1c/analysis/2015/qg_model/Jan17_drag_1e-3/vor_upperlayer.mp4', writer=writer, dpi=600)