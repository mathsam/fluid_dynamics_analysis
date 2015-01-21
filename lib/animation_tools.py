import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def make_animation(frames, num_frames):
    """
    make an animation with `num_frames` number of frames
    @param frames AnimationFrames object
    @param num_frames total number of frames one want to show
    @return fig, ani
    """
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    def animate_func(i):
        im = ax.imshow(frames.get_frame(i), cmap='gray')
        return im

    ani = animation.FuncAnimation(fig, animate_func, num_frames, interval=30)    
    return fig, ani
    
class AnimationFrames(object):
    def __init__(self, ncchain_raw, z_level):
        """
        asumme the input to be qg_model output
        @param ncchain_raw NetCDFChain object; or numpy array
        @param z_level which vertical level to output
        """
        self._ncchain = ncchain_raw
        self.total_frames = ncchain_raw.total_time_steps
        self._z_level = z_level
        return
    def get_frame(self, frame_index):
        import qg_transform
        psi_raw = self._ncchain[frame_index,:,:,:,self._z_level]
        psi_raw.shape += (1,)
        psic    = qg_transform.real2complex(psi_raw)
        vorc    = qg_transform.get_vorticity(psic)
        vorg    = qg_transform.spec2grid(vorc)
        vorg.shape = vorg.shape[0:2]
        return vorg