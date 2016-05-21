'''
Created on May 15, 2016

@author: Werner
'''
from element import Element
import numpy as np

def get_triangle_params(x, y):
        rot = lambda l, n: l[n:] + l[:n]
        ids = [0, 1, 2]
        a = [x[j] * y[k] - x[k] * y[j] for i, j, k in [rot(ids, n) for n in ids]]
        b = [y[j] - y[k] for i, j, k in [rot(ids, n) for n in ids]]
        c = [x[k] - x[j] for i, j, k in [rot(ids, n) for n in ids]]
        A = 0.5 * np.linalg.det(np.transpose(np.vstack((np.ones(3), x, y))))
        #print 'Triangle centre of gravity: x = %g, y = %g' %(cog_x,cog_y)
        #for i in range(3):
        #    print 'L%d = %g' %(i+1,Lxy(i,cog_x,cog_y))

        #bT = np.array([[i] for i in b])
        #cT = np.array([[i] for i in c])
        return (a,b,c, A)
    
class ElemPlateTriangleLin(Element):
    stiffdef = 'linear'
    stifftype_dict = {'linear':'stiffcalc_linear'}
    load = 0.0
    
    def get_load(self):
        p = np.zeros(9)
        #for i in range(3): p[i*3] = self.det_J[i] * self.load
        return p
    
    def plot(self, ax, magnitude, id=0):
        dom = self.domain()
        nds = self.nodes + [self.nodes[0]]
        xy = np.transpose(np.array([[dom.nodes[i - 1].x, dom.nodes[i - 1].y] for i in nds]))
        w = magnitude * np.transpose(np.hstack([dom.nodes[i - 1].v_disp for i in nds]))
        w0 = np.zeros_like(w)

        # Plot reference shape
        ax.plot_trisurf(xy[0, :], xy[1, :],  w[[i+id for i in range(0,12,3)]], color='green') #cmap=cm.jet)
        
    def stiffcalc_linear(self):
        nu = self.nu
        E = self.E
        nds = self.nodes

        # Bending thickness
        Cb = E*self.thickness**3/(1-nu**2)/12.*np.array([[1, nu, 0.],
                  [nu, 1, 0.],
                  [0., 0., (1-nu)/2.]])
        # Shear thickness
        Cs = 5*E*self.thickness/12./(1+nu) * np.eye(2)
        
        self.Cb = Cb
        self.Cs = Cs
        x=[self.domain().nodes[i - 1].x for i in nds]
        y=[self.domain().nodes[i - 1].y for i in nds]
        a,b,c,A = get_triangle_params(x,y)
        
        # Coordinates of the centre of gravity
        cgx = sum(x)/3.
        cgy = sum(y)/3.
        L = lambda i, x=cgx, y=cgy: (a[i] + b[i]*x + c[i]*y)/2./A
        
        Bb = np.zeros((3,9))
        for i in [2,5,8]:
            id = (i-2)/3
            Bb[0,i] = b[id]
            Bb[2,i-1] = -b[id]
            Bb[1,i-1] = -c[id]
            Bb[2,i] = c[id]
        Bb *= (1./2./A)
        Bs = np.zeros((2,9))
        for i in [2,5,8]:
            id = (i-2)/3
            Bs[0,i] = L(id)
            Bs[1,i-1] = -L(id)
            Bs[0,i-2] = b[id] / 2. / A
            Bs[1,i-2] = c[id] / 2. / A
            
        Kb = np.dot(np.transpose(Bb),np.dot(Cb,Bb)) * A
        Ks = np.dot(np.transpose(Bs),np.dot(Cs,Bs)) * A
        self.K = Kb + Ks

#get_triangle_params([0.,1.,5.],[0.,0.,3.])