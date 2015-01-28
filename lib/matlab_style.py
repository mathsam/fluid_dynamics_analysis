import numpy as np

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