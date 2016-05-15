'''
Created on May 15, 2016

@author: Werner
'''

from element import Element
import numpy as np

def get_triangle_stiffness_matrix(c1,c2, x, y):
        rot = lambda l, n: l[n:] + l[:n]
        ids = [0, 1, 2]
        a = [x[j] * y[k] - x[k] * y[j] for i, j, k in [rot(ids, n) for n in ids]]
        b = [y[j] - y[k] for i, j, k in [rot(ids, n) for n in ids]]
        c = [x[k] - x[j] for i, j, k in [rot(ids, n) for n in ids]]
        A = 0.5 * np.linalg.det(np.transpose(np.vstack((np.ones(3), x, y))))

        K_c1 = c1*A/12. * np.array([[(i==j)+1 for j in range(3)] for i in range(3)])
        bT = np.array([[i] for i in b])
        cT = np.array([[i] for i in c])
        K_c2 = c2/4./A * (bT*b + cT*c)
        K = K_c1 + K_c2
        return (K, A)
    
class ElemSoilQuadrangleLin(Element):
    stiffdef = '2triangles'
    stifftype_dict = {'2triangles':'set_K_2_triangles', '4triangles':'set_K_4_triangles',
                                   'isoparametric':'set_K_isoparametric'}

    def set_codes(self):
        '''
        Only the deflection dofs are considered
        '''
        self.v_code = [self.domain().nodes[i-1].v_code[0] for i in self.nodes]

    def set_K_2_triangles(self):
        nds = self.nodes + [self.nodes[0]]
        x = np.array([self.domain().nodes[i - 1].x for i in nds])
        y = np.array([self.domain().nodes[i - 1].y for i in nds])

        data = []
        data.append(get_triangle_stiffness_matrix(self.c1,self.c2, x[0:3], y[0:3]))
        data.append(get_triangle_stiffness_matrix(self.c1,self.c2, x[[2, 3, 0]], y[[2, 3, 0]]))
        k_lst, A_lst = zip(*data)

        K = np.zeros((4, 4))
        K[:3, :3] += k_lst[0]

        K[2:, 2:] += k_lst[1][:2, :2]
        K[:1, :1] += k_lst[1][2:, 2:]
        K[:1, 2:] += k_lst[1][2:, :2]
        K[2:, :1] += k_lst[1][:2, 2:]
        self.K = K
        self.A = sum(A_lst)

    def set_K_4_triangles(self):

        nds = self.nodes
        x = np.array([self.domain().nodes[i - 1].x for i in nds])
        y = np.array([self.domain().nodes[i - 1].y for i in nds])

        # Auxiliary center node
        xs = .25 * sum(x)
        ys = .25 * sum(y)

        # Triangle nodes
        xt_lst = [np.hstack((x[[i, (i + 1) % 4]], xs)) for i in range(4)]
        yt_lst = [np.hstack((y[[i, (i + 1) % 4]], xs)) for i in range(4)]

        data = [get_triangle_stiffness_matrix(self.c1,self.c2, xt, yt) for xt, yt in zip(xt_lst, yt_lst)]
        k_lst, A_lst = zip(*data)

        K = np.zeros((5, 5))
        for i in range(4):
            ids = [i, (i + 1) % 4, 4]
            for jj, j in enumerate(ids):
                for kk, k in enumerate(ids):
                    K[j, k] += k_lst[i][jj, kk]
        for i in range(4):
            K[i, :] -= K[i, 4] / K[4, 4] * K[4, :]
        self.K = K[:4, :4]
        self.A = sum(A_lst)

    def set_K_isoparametric(self):
        pass

    def plot(self, ax, magnitude, id=0):
        dom = self.domain()
        nds = self.nodes + [self.nodes[0]]
        xy = np.transpose(np.array([[dom.nodes[i - 1].x, dom.nodes[i - 1].y] for i in nds]))
        w = magnitude * np.transpose(np.array([dom.nodes[i - 1].v_disp[0] for i in nds]))
        w0 = np.zeros_like(w)

        # Plot reference shape
        ax.plot_trisurf(xy[0, :], xy[1, :],  w[[i+id for i in range(len(w))]], color='blue') #cmap=cm.jet)