'''
Created on May 12, 2016

@author: Werner
'''

from FriendlyFEM import *
nelem = 10
length = 3.
load = [0., 1., 0.]
nodes = []
elements = []

for i in range(nelem + 1):
    nodes.append(Node(x=i * length / float(nelem), y=0., F_ext=[0., 0., 0.], supports=[0, 0, 0], nnodedofs=3))
    if i < nelem:
        elements.append(ElemBeamLin(nodes=[i + 1, i + 2], A=.16, Iy=.0021, Ik=.0021, k=1., E=20000., nu=0.2, density=0.0))
nodes[0].supports = [1, 1, 1]
nodes[-1].F_ext = load

dom = Domain(nodes=nodes, elements=elements, c1=1., c2=1.)
dom.solve(verbose=True)
dom.plot(magnitude=1.)
