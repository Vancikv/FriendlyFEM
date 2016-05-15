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

    '''
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
    '''
                
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

base_optional = ['nodes', 'stiffdef']
base_compulsory = ['E', 'nu', 'density']
beam_lin_optional = base_optional + []
beam_lin_compulsory = base_compulsory + ['A', 'Iy', 'Ik', 'k']

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

class ElemBeamLin(Element):
    def __init__(self, *args, **kwargs):
        for key in beam_lin_compulsory:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                print 'Error: The property %s for a beam element has not been specified.' % key
                raise
        for key in beam_lin_optional:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                pass
            
    def plot(self, ax, magnitude, id=0):
        dom = self.domain()
        nds = self.nodes + [self.nodes[0]]
        xy = np.transpose(np.array([[dom.nodes[i - 1].x, dom.nodes[i - 1].y] for i in nds]))
        w = magnitude * np.transpose(np.array([dom.nodes[i - 1].v_disp[0] for i in nds]))
        w0 = np.zeros_like(w)

        # Plot reference shape
        ax.plot(xy[0, :], xy[1, :],  w0, color = 'black', linewidth=1.) #linestyle='-', color='black', linewidth=1.
        # Plot deformed shape
        ax.plot(xy[0, :], xy[1, :],  w, color = 'red', linewidth=5.) #linestyle='-', color='black', linewidth=1.
            
    def set_matrices(self):
        '''
        Approximating the deflection and two rotations (bending, torsion) with linear functions
        '''
        # Shear and bending part
        E, k, A, I, Ik = self.E, self.A, self.k, self.Iy, self.Ik
        G = E/2./(1+self.nu)
        n1 = self.domain().nodes[self.nodes[0]-1]
        n2 = self.domain().nodes[self.nodes[1]-1]
        L = ((n1.x-n2.x)**2+(n1.y-n2.y)**2)**0.5
        K_b = np.zeros((4,4))
        K_b[1,1] = K_b[3,3] = 1.
        K_b[3,1] = K_b[1,3] = -1.
        K_b *= E*I/L
        K_s = k*G*A* np.array([[1./L, -1./2., -1./L, -1./2.],
                               [-1./2., L/4., 1./2., L/4.],
                               [-1./L, 1./2., 1./L, 1./2.],
                               [-1./2., L/4., 1./2., L/4.],
                               ])
        # Torsion part
        K_t = G*Ik/L*np.array([[1.,-1.],
                               [-1.,1.]
                               ])
        K = np.zeros((6,6))
        K[:4,:4] += K_b + K_s
        K[4:,4:] += K_t
        # Swap rows to match vector (w1,fiy1,fix1,w2,fiy2,fix2)
        K[[2,3,4],:] = K[[4,2,3],:]
        K[:,[2,3,4]] = K[:,[4,2,3]]
        self.K_loc = K
        # Obtain transformation matrix
        s = (n2.y-n1.y) / L
        c = (n2.x-n1.x) / L
        self.T = np.array([[1.,0.,0.,0.,0.,0.],
                           [0.,c,-s,0.,0.,0.],
                           [0.,s,c,0.,0.,0.],
                           [0.,0.,0.,1.,0.,0.],
                           [0.,0.,0.,0.,c,-s],
                           [0.,0.,0.,0.,s,c],
                           ])
        self.K = np.dot(np.dot(np.transpose(self.T), K), self.T)
        
node_compulsory = ['x','y','nnodedofs','F_ext','supports']

class Node(object):
    def __init__(self, *args, **kwargs):
        for key in node_compulsory:
            try:
                setattr(self,key,kwargs[key])
            except KeyError:
                print 'Error: The property %s for a node element has not been specified.' % key
                raise
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
        
if __name__ == '__main__':
    nelem = 3
    length = 3.
    load = [-1.,0.,0.]
    nodes = []
    elements = []
    
    for i in range(nelem+1):
        nodes.append(Node(x=i*length/float(nelem), y=0., F_ext=[0.,0.,0.], supports=[0,0,0], nnodedofs=3))
        if i<nelem:
            elements.append(ElemBeamLin(nodes=[i+1,i+2],A=.16,Iy=.0021,Ik=.0021,k=1.,E=20000.,nu=0.2, density=0.0))
    nodes[0].supports = [1,1,1]
    nodes[-1].F_ext = load
    
    dom = Domain(nodes=nodes,elements=elements)
    dom.solve(verbose = True)
    dom.plot(magnitude=1.)