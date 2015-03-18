import numpy as np
import os
import matplotlib.pyplot as plt

def clear():
    """
    Emulate Matlab command clear
    
    Clean objects binded to __main__
    """
    import sys
    current = sys.modules['__main__']
    for var in dir(current):
        if var[0] != '_':
            current.__dict__.pop(var)
        
def save(filename, **kwargs):
    """
    Usage:
        save('mydata',varname1=var1, varname2=var2)
    """
    np.savez(filename,**kwargs)
    
def load(filename):
    """
    Load the numpy arrays in `filename` and bind them to __main__ namespace
    Usage:
        load('mydata')
    """
    import sys
    main_dict = sys.modules['__main__'].__dict__
    if len(filename)<5 or filename[-3:] != '.npz':
        filename += '.npz'
    var = np.load(filename)
    for key, value in var.iteritems():
        print '%s added to __main__ namespace' %key
        main_dict[key] = value
        
def fig_save(fig, filepath):
    """
    save a figure `fig` to filename and dir specified in `filepath`
    If directory does not exist, create it.
    @param fig figure object
    @param filepath absolute path or just a filename
    Example:
        fig_save(fig, '~/analysis/myexp/fig.png')
        fig_save(fig, 'fig.png') #save to current working dir
    """
    try:
        fig.savefig(filepath)
    except IOError:
        endindex = filepath[::-1].index('/')
        savedir = filepath[:endindex]
        if not os.path.isdir(savedir):
            os.mkdir(savedir)
        fig.savefig(filepath)
    return None
    
def plotyy(x1, y1, x2, y2, xlabel=None, ylabel1=None, ylabel2=None):
    """
    Create a figure with two y-axes
    """
    fig, ax1 = plt.subplots()
    ax1.plot(x1, y1, 'b-')
    if xlabel:
        ax1.set_xlabel(xlabel)
    if ylabel1:
        ax1.set_ylabel(ylabel1)
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
        
    ax2 = ax1.twinx()
    ax2.plot(x2, y2, 'r-')
    if ylabel2:
        ax2.set_ylabel(ylabel2)
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    return fig