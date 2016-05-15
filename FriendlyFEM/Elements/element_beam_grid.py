'''
Created on May 12, 2016

@author: Werner
'''

import numpy as np
from element import Element

class ElemBeamGrid(Element):
    stiffdef = 'linear'
    stifftype_dict = {'linear':'stiffcalc_linear'}
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
            
    def stiffcalc_linear(self):
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