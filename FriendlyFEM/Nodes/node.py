'''
Created on May 12, 2016

@author: Werner
'''
import numpy as np
from FriendlyFEM.Auxiliary import CommonInit

class Node(CommonInit):
    def __init__(self, **kwargs):
        super(Node,self).__init__(**kwargs)
        self.free_dofs = [i for i, j in enumerate(self.supports) if j == 0]
        self.v_disp = np.zeros(self.nnodedofs)
        
    def set_codes(self, maxcode, verbose=False):
        self.v_code = []
        for s in self.supports:
            if s == 0:
                maxcode += 1
                self.v_code.append(maxcode)
            else:
                self.v_code.append(0)
        if verbose: print self.v_code
        return maxcode

    def set_disp(self, val):
        self.v_disp[self.free_dofs] = val[self.free_dofs]