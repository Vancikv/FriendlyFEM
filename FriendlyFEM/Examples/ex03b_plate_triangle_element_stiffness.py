'''
Created on May 18, 2016

@author: Werner
'''

from FriendlyFEM import *

lx = 1.
ly = 1.
E = 20000.
nu = 0.2
thickness = 0.2
nodes = [Node(x=0., y=0., F_ext=[0., 0., 0.], supports=[1, 1, 1], nnodedofs=3),
         Node(x=lx, y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=2.*lx, y=0., F_ext=[-1., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=0., y=ly, F_ext=[0., 0., 0.], supports=[1, 1, 1], nnodedofs=3),
         Node(x=lx, y=ly, F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
        Node(x=2.*lx, y=ly, F_ext=[-1., 0., 0.], supports=[0, 0, 0], nnodedofs=3)]

elements = [ElemPlateTriangleLin(nodes=[1,2,5], E=E, nu=nu, density=0.0, thickness=thickness),
            ElemPlateTriangleLin(nodes=[1,5,4], E=E, nu=nu, density=0.0, thickness=thickness),
            ElemPlateTriangleLin(nodes=[2,3,6], E=E, nu=nu, density=0.0, thickness=thickness),
            ElemPlateTriangleLin(nodes=[2,6,5], E=E, nu=nu, density=0.0, thickness=thickness)]

dom = Domain(elements=elements, nodes = nodes)
dom.solve(verbose=True)
K = dom.K_glob
print np.linalg.eigvals(K)