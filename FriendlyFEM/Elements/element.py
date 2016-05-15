'''
Created on May 12, 2016

@author: Werner
'''
import numpy as np
import itertools as it
from FriendlyFEM.Auxiliary import CommonInit

class Element(CommonInit):
    domain = None
    nodes = []
    stifftype_dict = {}
    stiffdef = ''

    def calculate_stiffness(self):
        '''
        Calculate and store the element stiffness matrix. 
        The stiffness definition is specified by the variable stiffdef.
        '''
        getattr(self, self.stifftype_dict[self.stiffdef])()
    
    def get_load(self):
        return np.zeros_like(self.v_code)
        
    def set_codes(self):
        self.v_code = list(it.chain(*[self.domain().nodes[i-1].v_code for i in self.nodes]))

    def plot(self, ax, magnitude):
        pass
