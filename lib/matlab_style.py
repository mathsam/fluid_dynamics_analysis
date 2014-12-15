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
