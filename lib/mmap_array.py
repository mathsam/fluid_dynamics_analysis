import numpy as np

class MmapArray(np.ndarray):
    """
    work the sampe as numpy.ndarray but has an additional attribute, which is a
    file handler the array is memory mapped to
    """
    
    def __new__(cls, input_array, fh):
        """
        Initialized an array with a file handler.
        
        Args:
            input_array: numpy.ndarray
            fh: file handler the array is memory mapped to
        """
        obj = np.asarray(input_array).view(cls)
        obj.fh = fh
        return obj
        
    def __array_finalize__(self, obj):
        if obj is None: return
        self.obj = getattr(obj, 'fh', None)