############################################################################################
# External includes

import numpy
from collections.abc import MutableMapping as MM

############################################################################################
# Custom dictionary class

class fixedSizeDict ( MM ) :

    _internal = []
    
    def __init__ ( self, Nobj, *args, **kwargs ) :

        self.size = Nobj

        for arg in args :
            if isinstance( arg, MM ) :
                for k,v in arg.items() :
                    self[k] = v
                
        for k, v in kwargs.items() :
            self[k] = v
            
    def __getitem__ ( self, key ) :
        return self.__dict__[ key ]
        
    def __setitem__ ( self, key, value ) :

        if hasattr( value, '__len__' ) :
            if len( value ) == self.size :
                try :
                    self.__dict__[key] = numpy.array( value )
                    # setattr( self, key,
                    #          numpy.array( value ) )
                except :
                    raise
            else :
                raise TypeError( f'All fields should have the same lenght '
                                 f'(={self.size}) but '
                                 f'len({key})={len(value)}!={self.size}' )
        else :
            raise TypeError( f"'{key}' value not valid: {type(self).__name__} "
                             f"values should be iterables of lenght={self.size}" )
        
    def __delitem__ ( self ) :

        raise RuntimeError( "Deletion not allowed" ) 

    def __iter__ ( self ) :
        return iter( self.__dict__ )
            
    def size ( self ) :
        return self.size
    
    def __len__ ( self ) :
        return self.size
    
    def __repr__ ( self ) :
        return f"{type(self).__name__}({repr( self.__dict__ )})"
    
    def __str__ ( self ) :
        _str = f"{type(self).__name__}( "
        pad = " " * len(_str)
        for k,v in self.__dict__.items() :
            if k not in type(self)._internal :
                _str += f"\n{pad}{k} = {str(v.data)}"
        _str += "\n{pad})"
        return _str
    
    def __getstate__ ( self ) :
        """ Allow pickle
        """
        return self.__dict__
    
    def __setstate__ ( self, d ) :
        """ Allow unpickle
        """
        self.__dict__ = d 
        for k, v in d.items() :
            setattr( self, k, v )
        return;
    
    def rename_field ( self, old_key, new_key ) :
        
        if old_key in self.__dict__.keys() :
            self.__setitem__( new_key, self.__dict__.pop( old_key ) )
            return;
        
    def keys ( self ) :
        return self.__dict__.keys()
    
    def values ( self ) :
        return self.__dict__.values()
    
    def items ( self ) :
        return self.__dict__.items()
        
    def pop ( self, key ) :
        if key in self.__dict__.keys() :
            value = self.__dict__.pop( key )
            return value
        return None

    def has_key ( self, key ) :
        if key not in self.keys() :
            return False
        return True
    
############################################################################################

class maskedDict ( MM ) :

    _internal = [ 'mask' ]
    
    def __init__ ( self, Nobj, *args, **kwargs ) :
        import numpy as np

        self.mask = np.zeros( ( Nobj, ), dtype = bool )

        for arg in args :
            if isinstance( arg, MM ) :
                for k,v in arg.items() :
                    self.__setitem__( k, v )
                
        for k, v in kwargs.items() :
            self.__setitem__( k, v )
            
    def __getitem__ ( self, key ) :
        return self.__dict__[ key ].data
        
    def __setitem__ ( self, key, value ) :
        import numpy as np

        if key == 'mask' : return 
        if hasattr( value, '__len__' ) :
            if len( value ) == len( self.mask ) :
                try :
                    setattr( self, key,
                             np.ma.masked_array( value, mask = self.mask, hard_mask = True ) )
                except np.ma.MaskError as merr :
                    raise np.ma.MaskError( f"Setting key '{key}' raised MaskError:\n{merr}" )
            else :
                raise TypeError( f'Masked fields should all have the same lenght '
                                 f'(={len(self.mask)}) but '
                                 f'len({key})={len(value)}!={len(self.mask)}' )
        else :
            raise TypeError( f"'{key}' value not valid: {type(self).__name__} "
                             f"values should be iterables of lenght={len(self.mask)}" )
        
    def __delitem__ ( self ) :

        raise RuntimeError( "Deletion not allowed" ) 

    def __iter__ ( self ) :
        return iter( self.__dict__ )
            
    def size ( self ) :
        return len(self.mask)
    
    def __len__ ( self ) :
        return self.size() - self.mask.sum()
    
    def __repr__ ( self ) :
        return f"{type(self).__name__}({repr( self.__dict__ )})"
    
    def __str__ ( self ) :
        _str = f"{type(self).__name__}( "
        pad = " " * len(_str)
        _str += f"mask = {str(self.mask)}"
        for k,v in self.__dict__.items() :
            if k not in type(self)._internal :
                _str += f"\n{pad}{k} = {str(v.data)}"
        _str += " )"
        return _str
    
    def __getstate__ ( self ) :
        """ Allow pickle
        """
        return self.__dict__
    
    def __setstate__ ( self, d ) :
        """ Allow unpickle
        """
        self.__dict__ = d 
        for k, v in d.items() :
            setattr( self, k, v )
        return;
    
    def rename_field ( self, old_key, new_key ) :
        
        if old_key in self.__dict__.keys() :
            self.__setitem__( new_key, self.__dict__.pop( old_key ) )
            return;
        
    def keys ( self ) :
        return self.__dict__.keys()
    
    def values ( self ) :
        return self.__dict__.values()
    
    def items ( self ) :
        return self.__dict__.items()
        
    def pop ( self, key ) :
        if key in self.__dict__.keys() :
            value = self.__dict__.pop( key )
            return value
        return None

    def has_key ( self, key ) :
        if key not in self.keys() :
            return False
        return True

    def reset_mask ( self ) :
        import numpy as np
        self.mask *= np.zeros_like( self.mask )
        return;

    def add_mask ( self, mask ) :
        self.mask += mask
        return;

    def iter_masked ( self, *keys ) :

        for ii in zip( *[ self.__dict__[k][~self.mask] for k in keys ] ) :
            yield ii
    
############################################################################################
