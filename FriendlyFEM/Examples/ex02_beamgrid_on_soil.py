'''
Created on May 15, 2016

@author: Werner
'''

from FriendlyFEM import *

dom = beamgrid_on_soil(N_lst=[[0.,0.],[10.,0.],[10.,10.],[0.,10.]],
               n_beam_1=3, n_div_1=4, n_out_1=3, d_out_1=0.5, 
               n_beam_2=3, n_div_2=4, n_out_2=3, d_out_2=0.5,
               A=.16,Iy=0.0021,Ik=.0021,k=1.,E=20000., c1=10000.,c2=1.e6)
def find_node_by_coords(x,y,domain):
    isclose = lambda n1, n2, tol = 1e-3: (((n2[0] - n1[0]) ** 2 + (n2[1] - n1[1]) ** 2) ** 0.5) < tol
    for i, nd in enumerate(domain.nodes):
        if isclose((nd.x,nd.y),(x,y)):
            return i
for x,y in [[0.,0.],[10.,0.],[10.,10.],[0.,10.],[5.,5.]
            #[10./3.,10./3.],[10./3.,20./3.],[20./3.,10./3.],[20./3.,20./3.],
            ]:
    n = find_node_by_coords(x,y,dom)
    dom.nodes[n].F_ext = [-500.,0.,0.]


dom.solve()
dom.plot()
