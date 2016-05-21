'''
Created on May 18, 2016

@author: Werner
'''

from FriendlyFEM import *

nodes = [Node(x=0., y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=1., y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=1., y=1., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3)]

elements = [ElemPlateTriangleLin(nodes=[1,2,3], E=20000., nu=0.2, density=0.0, thickness=0.2)]

dom = Domain(elements=elements, nodes = nodes)

print 'Single triangle stiffness matrix eigenvalues:'
dom.global_assembly()
K = dom.K_glob
print np.linalg.eigvals(K)

nodes = [Node(x=0., y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=1., y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
         Node(x=1., y=1., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3),
        Node(x=0., y=1., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3)]

elements = [ElemPlateTriangleLin(nodes=[1,2,3], E=20000., nu=0.2, density=0.0, thickness=0.2),
            ElemPlateTriangleLin(nodes=[1,3,4], E=20000., nu=0.2, density=0.0, thickness=0.2)]

dom = Domain(elements=elements, nodes = nodes)

print '\nTwo connected triangles\' stiffness matrix eigenvalues:'
dom.global_assembly()
K = dom.K_glob
print np.linalg.eigvals(K)