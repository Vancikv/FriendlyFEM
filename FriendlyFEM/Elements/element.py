'''
Created on May 12, 2016

@author: Werner
'''
import numpy as np

class Element(object):
    domain = None
    nodes = []

    def set_matrices(self):
        '''
        Calculate and store local stiffness, mass and damping matrices
        '''
        pass
    
    def get_load(self):
        return np.zeros_like(self.v_code)
        
    def set_codes(self):
        self.v_code = reduce(lambda x,y: x+y, [self.domain().nodes[i-1].v_code for i in self.nodes])

    def plot(self, ax, magnitude):
        pass
