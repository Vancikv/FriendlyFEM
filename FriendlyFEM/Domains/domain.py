'''
Created on May 12, 2016

@author: Werner
'''

import numpy as np
import weakref
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools as it
domain_optional = ['elements', 'nodes']

class Domain(object):
    @property
    def elements(self):
        return self._elements
    @elements.setter
    def elements(self, value):
        self._elements = value
        for e in self._elements:
            e.domain = weakref.ref(self)
        
    nodes = []
    def __init__(self, *args, **kwargs):
        for key in domain_optional:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                pass

    def plot(self, id=None,magnitude=10.):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for el in self.elements:
            if id is not None:
                el.plot(ax, magnitude,id=id)
            else:
                el.plot(ax, magnitude)
        #ax.set_zlim3d(-.2)
        plt.show()

    def node_elems(self, n):
        cnt = 0
        for el in self.elements:
            if (n+1) in el.nodes:
                cnt += 1
        return cnt
    
    def id_by_coords(self, x,y):
        isclose = lambda n1, n2, tol = 1e-6: (((n2[0] - n1[0]) ** 2 + (n2[1] - n1[1]) ** 2) ** 0.5) < tol
        for i, nd in enumerate(self.nodes):
            if isclose((nd.x,nd.y),(x,y)):
                return i
        return None
        
    def global_assembly(self):
        # Initialize local values
        els = self.elements
        nds = self.nodes
        code_count = 0
        for nd in nds:
            code_count = nd.set_codes(code_count)
        print "code count = %d" % code_count
        self.code_count = code_count

        K_glob = np.zeros((code_count, code_count))  # Global stiffness
        for el in els:
            el.set_matrices()
            el.set_codes()
            cds = [(i, el.v_code[i]) for i in range(len(el.v_code)) if el.v_code[i] != 0]
            K_loc = el.K
            for ii, ij in it.product(cds, cds):
                K_glob[ii[1] - 1, ij[1] - 1] += K_loc[ii[0], ij[0]]
        self.K_glob = K_glob    
            
        f_glob = np.zeros(code_count)  # Global load vector
        for nd in nds:
            for i, c in enumerate(nd.v_code):
                if c != 0: f_glob[c - 1] += nd.F_ext[i]
        for el in els:
            p = el.get_load()
            for i, c in enumerate(el.v_code):
                if c != 0: f_glob[c - 1] += p[i]
        self.f_glob = f_glob
                
    def solve(self, verbose=False):
        self.global_assembly()
        K_glob = self.K_glob
        f_glob = self.f_glob

        # Call numpy linear solver
        u_res = np.linalg.solve(K_glob,f_glob)
        
        if verbose:
            K_inv = np.linalg.inv(K_glob)
            print 'Global stiffness matrix:', K_glob
            print 'Inverse stiffness matrix:', K_inv
            print 'Control product K * K_inv:', np.round(np.dot(K_glob,K_inv),4)
            print "Global load vector", f_glob
            print "Resulting displacement vector", u_res
            
        # Distribute results
        for nd in self.nodes:
            nd.set_disp(u_res[[i-1 for i in nd.v_code]])
