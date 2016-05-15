'''
Created on May 15, 2016

@author: Werner
'''
from FriendlyFEM.Elements import ElemBeamGrid, ElemSoilQuadrangleLin
from FriendlyFEM.Nodes import Node
from FriendlyFEM.Domains import Domain
import numpy as np

def beamgrid_on_soil(N_lst,n_beam_1, n_div_1, n_out_1, d_out_1, n_beam_2, n_div_2, n_out_2, d_out_2,
               A=.16,Iy=.0021,Ik=.0021,k=1.,E=20000., density=0.,nu=0.2,c1=1000.,c2=1000.):

    if len(N_lst) != 4:
        print 'Wrong number of nodes for domain generation'
        return
    # Create side direction and normal vectors
    x = np.array(N_lst+[N_lst[0]])
    v = np.array([x[i+1,:]-x[i,:] for i in range(4)])
    n = np.array([[i[1] * 1./np.linalg.norm(i),-i[0] * 1./np.linalg.norm(i)] for i in v])
    
    # Step 1 - core:
    ndiv1 = (n_beam_1-1)*n_div_1+1 
    ndiv2 = (n_beam_2-1)*n_div_2+1
    
    get_edge_division = lambda N1, N2, ndiv, endp=True: np.transpose(np.vstack( (np.linspace(N1[0],N2[0],ndiv,endpoint=endp),np.linspace(N1[1],N2[1],ndiv,endpoint=endp)) ))
    
    col_edge_1 = get_edge_division(x[3],x[0], ndiv2)
    col_edge_2 = get_edge_division(x[2],x[1], ndiv2)
    N_core = np.array([get_edge_division(ce1,ce2,ndiv1) for ce1, ce2 in zip(col_edge_1, col_edge_2)])
    
    # Step 2 - bottom perimeter: vl=v4, vr=v2, vb=v1, N1=N1, N2=N2, nn=n1
    vl, vr, vb, N1, N2, nn = v[3], v[1], v[0], x[0], x[1], n[0]
    Nb = N1 + n_out_2 * d_out_2 * nn
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Matrix of line intersection equation
    rhs = Nb - N1 # Right-hand side ...
    res = np.dot(np.linalg.inv(Mbl),rhs) # Result - line equation parameters
    Nl = N1 + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) )) # The same for the other point
    rhs = Nb - N2
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nr = N2 + res[0] * vr
    
    edge_1 = get_edge_division(N1,Nl, n_out_2)
    edge_2 = get_edge_division(N2,Nr, n_out_2)
    N_bottom = np.array([get_edge_division(ce1,ce2,ndiv1) for ce1, ce2 in zip(edge_1, edge_2)])
    N_bottom = N_bottom[1:,:]
    
    # Save the new corner points - needed for left and right
    Nbl = Nl
    Nbr = Nr

    # Step 2 - top perimeter: vl=v2, vr=v4, vb=v3, N1=N3, N2=N4, nn=n3
    vl, vr, vb, N1, N2, nn = v[1], v[3], v[2], x[2], x[3], n[2]
    Nb = N1 + n_out_2 * d_out_2 * nn
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Matrix of line intersection equation
    rhs = Nb - N1 # Right-hand side ...
    res = np.dot(np.linalg.inv(Mbl),rhs) # Result - line equation parameters
    Nl = N1 + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) )) # The same for the other point
    rhs = Nb - N2
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nr = N2 + res[0] * vr
    
    edge_1 = get_edge_division(Nr, N2, n_out_2)
    edge_2 = get_edge_division(Nl, N1, n_out_2)
    N_top = np.array([get_edge_division(ce1,ce2,ndiv1) for ce1, ce2 in zip(edge_1, edge_2)])
    N_top = N_top[:-1,:]

    Ntl = Nr
    Ntr = Nl

    # Step 3 - left perimeter: vl=v3, vr=v1, vb=v4, N1=Ntl, N2=Nbl, nn=n4
    vl, vr, vb, N1, N2, nn = v[2], v[0], v[3], Ntl, Nbl, n[3]
    Nb = N1 + n_out_2 * d_out_2 * nn
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Matrix of line intersection equation
    rhs = Nb - N1 # Right-hand side ...
    res = np.dot(np.linalg.inv(Mbl),rhs) # Result - line equation parameters
    Nl = N1 + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) )) # The same for the other point
    rhs = Nb - N2
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nr = N2 + res[0] * vr
    
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Intermediate points
    rhs = Nb - x[3]
    res = np.dot(np.linalg.inv(Mbl),rhs)
    Nli = x[3] + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) ))
    rhs = Nb - x[0]
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nri = x[0] + res[0] * vr

    edge_1 = np.vstack( (get_edge_division(Nl,Nli, n_out_2-1, endp=False),
                         get_edge_division(Nli,Nri, ndiv2-1, endp=False),
                         get_edge_division(Nri,Nr, n_out_2)) )
    edge_2 = np.vstack( (get_edge_division(N1,x[3], n_out_2-1, endp=False),
                         get_edge_division(x[3],x[0], ndiv2-1, endp=False),
                         get_edge_division(x[0],N2, n_out_2)) )
    N_left = np.array([get_edge_division(ce1,ce2,n_out_1) for ce1, ce2 in zip(edge_1, edge_2)])
    N_left = N_left[:,:-1]

    # Step 3 - right perimeter: vl=v1, vr=v3, vb=v2, N1=Nbr, N2=Ntr, nn=n2
    vl, vr, vb, N1, N2, nn = v[0], v[2], v[1], Nbr, Ntr, n[1]
    Nb = N1 + n_out_2 * d_out_2 * nn
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Matrix of line intersection equation
    rhs = Nb - N1 # Right-hand side ...
    res = np.dot(np.linalg.inv(Mbl),rhs) # Result - line equation parameters
    Nl = N1 + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) )) # The same for the other point
    rhs = Nb - N2
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nr = N2 + res[0] * vr
    
    Mbl = np.transpose(np.vstack( (vl,-vb) )) # Intermediate points
    rhs = Nb - x[1]
    res = np.dot(np.linalg.inv(Mbl),rhs)
    Nli = x[1] + res[0] * vl
    
    Mbr = np.transpose(np.vstack( (vr,-vb) ))
    rhs = Nb - x[2]
    res = np.dot(np.linalg.inv(Mbr),rhs)
    Nri = x[2] + res[0] * vr

    edge_1 = np.vstack( (get_edge_division(N2,x[2], n_out_2-1, endp=False),
                         get_edge_division(x[2],x[1], ndiv2-1, endp=False),
                         get_edge_division(x[1],N1, n_out_2)) )
    edge_2 = np.vstack( (get_edge_division(Nr,Nri, n_out_2-1, endp=False),
                         get_edge_division(Nri,Nli, ndiv2-1, endp=False),
                         get_edge_division(Nli,Nl, n_out_2)) )
    N_right = np.array([get_edge_division(ce1,ce2,n_out_1) for ce1, ce2 in zip(edge_1, edge_2)])
    N_right = N_right[:,1:]
    
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(N_core[:,:,0], N_core[:,:,1], marker='o', color='black')
    ax.plot(N_bottom[:,:,0], N_bottom[:,:,1], marker='o', color='green')
    ax.plot(N_top[:,:,0], N_top[:,:,1], marker='o', color='green')
    ax.plot(N_left[:,:,0], N_left[:,:,1], marker='o', color='blue')
    ax.plot(N_right[:,:,0], N_right[:,:,1], marker='o', color='blue')
    plt.axis('equal')
    plt.show()
    '''
    coords = np.hstack( (N_left, 
                         np.vstack( (N_top, N_core, N_bottom) ),
                         N_right) )
    nodes = []
    elements = []
    dimJ = coords.shape[1]
    dimI = coords.shape[0]
    node_n = lambda i,j: i*dimJ + j + 1
    elem_nodes = lambda i,j: [node_n(i+1,j),node_n(i+1,j+1),node_n(i,j+1),node_n(i,j)]
    beam_grid_i = [n_out_2-1 + i*n_div_2 for i in range(n_beam_2)]
    beam_grid_j = [n_out_1-1 + i*n_div_1 for i in range(n_beam_1)]
     
    for i in range(dimI):
        for j in range(dimJ):
            if ((i in beam_grid_i) and (beam_grid_j[0]<j<beam_grid_j[-1])):
                nodes.append(Node(x=coords[i,j,0],y=coords[i,j,1], nnodedofs=3, F_ext=[0.,0.,0.], supports=[0,0,0]))
            else:
                nodes.append(Node(x=coords[i,j,0],y=coords[i,j,1], nnodedofs=1, F_ext=[0.], supports=[0]))
            if j<(dimJ-1) and i<(dimI-1):
                elements.append(ElemSoilQuadrangleLin(nodes=elem_nodes(i,j),c1=c1,c2=c2,E=0.,nu=0.,density=0.))
                
    for i in beam_grid_i:
        for j in range(beam_grid_j[0],beam_grid_j[-1]):
            elements.append(ElemBeamGrid(nodes=[node_n(i,j),node_n(i,j+1)],A=A,Iy=Iy,Ik=Ik,k=k,E=E, density=density,nu=nu))
    for j in beam_grid_j:
        for i in range(beam_grid_i[0],beam_grid_i[-1]):
            elements.append(ElemBeamGrid(nodes=[node_n(i,j),node_n(i+1,j)],A=A,Iy=Iy,Ik=Ik,k=k,E=E, density=density,nu=nu))
    return Domain(nodes=nodes,elements=elements)
