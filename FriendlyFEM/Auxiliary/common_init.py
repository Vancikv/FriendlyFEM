'''
Created on May 15, 2016

@author: Werner

The base class defined here defines the unified initialization style
'''
from argument_library import *

class CommonInit(object):
    #def __new__(cls,*args, **kwargs):
    #    return super(CommonInit,cls).__new__(cls,*args, **kwargs)
        
    def __init__(self, **kwargs):
        nm = self.__class__.__name__
        cmp_args = globals()['cmp_'+nm]
        opt_args = globals()['opt_'+nm]
        for key in cmp_args:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                raise RuntimeError('The property %s has not been specified for object of type %s' % (key, self.__class__))
                
        for key in opt_args:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                pass