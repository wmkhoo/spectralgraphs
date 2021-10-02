import matplotlib.pyplot as plt
import networkx as nx
import scipy
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
import numpy as np


def house_graphs():
    g = nx.Graph()
    g.add_edges_from([(1,2),(1,5),(1,7),(2,6),(3,4),(3,5),(3,7),(4,6),(5,6),(5,7),(5,8)])
    h = nx.Graph()
    h.add_edges_from([(1,2),(1,5),(2,6),(2,7),(3,4),(3,5),(4,6),(4,7),(5,6),(6,7),(6,8)])
    return g,h


# Finds diagonal matrix S, complexity 2^n
# Input: graphs g, h
# Output: permutation matrix "P", minimised distance "min"
def min_P(g, h):
    A = nx.adjacency_matrix(g)
    B = nx.adjacency_matrix(h)
    wA, vA = eigh(A.todense())
    wB, vB = eigh(B.todense())
    min = 999
    P_min = None
    n = len(A.todense())
    for i in range(2 ** n):
        m = i
        a = []
        for j in range(n):
            if m & 1:
                a.append(1)
            else:
                a.append(-1)
            m = m >> 1
        S = np.diag(a)
        P = vA @ S @ vB.T
        norm = np.linalg.norm(B - (P.T @ A @ P))
        if norm < min:
            P_min = P
            min = norm
    return P_min, min


# Algorithm by Stefan Klus and Tuhin Sahai, "A Spectral Assignment Approach for the Graph Isomorphism Problem", https://arxiv.org/pdf/1411.0969.pdf
# Input: Networkx graphs g & h
# Output: Mapping of nodes in g to h
# TODO: Doesn't yet deal with repeated eigenvalues
def min_P_la(g, h):
    # Make sure g and h have the same number of nodes
    if len(g.nodes()) != len(h.nodes()):
        return []
    
    A = nx.adjacency_matrix(g)
    B = nx.adjacency_matrix(h)
    wA, vA = eigh(A.todense())
    wB, vB = eigh(B.todense())
    n = len(A.todense())
    
    for i in range(n):
        # For ambiguous vectors (elements sum to 0), take the absolute values
        if abs(vA[:,i].sum()) < 0.00000001:
            vA[:,i] = abs(vA[:,i])
        # Make sure 1.T v[i] > 0
        elif vA[:,i].sum() < 0:
            vA[:,i] = -1 * vA[:,i]
        
        # For ambiguous vectors (elements sum to 0), take the absolute values
        if abs(vB[:,i].sum()) < 0.00000001:
            vB[:,i] = abs(vB[:,i])
        # Make sure 1.T v[i] > 0
        elif vB[:,i].sum() < 0:
            vB[:,i] = -1 * vB[:,i]
    
    C = []
    for i in range(n):
        c_row = []
        for j in range(n):
            # c_ij = Sum_k | v[i] - w[i] |
            c_row.append(abs((vA[i] - vB[j]).sum()))
        C.append(c_row)
    C = np.array(C)
    
    # c = min (P.T @ C), i.e. P is the solution to linear assignment of C
    row_ind, col_ind = linear_sum_assignment(C)
    
    # Compute the node mapping from g to h
    mapping = []
    ns1 = list(g.nodes())
    ns2 = list(h.nodes())
    for i in range(len(row_ind)):
        mapping.append((ns1[row_ind[i]], ns2[col_ind[i]]))
    return mapping


def plot_eigen(x):
    plt.figure(figsize=(8,6))
    plt.scatter(np.arange(len(x)), x, marker="+")
    plt.title("Eigenvalues plot")
    plt.show()


def plot_eigen2(x1, x2):
    plt.figure(figsize=(8,6))
    plt.scatter(np.arange(len(x1)), x1, marker="x")
    plt.scatter(np.arange(len(x2)), x2, marker="+")
    plt.title("Eigenvalues plot")
    plt.show()


def label2color(labels):
    colormap = []
    for l in labels:
        if l == 0:
            colormap.append('blue')
        elif l == 1:
            colormap.append('green')  
        elif l == 2:
            colormap.append('red') 
        elif l == 3: 
            colormap.append('yellow')
        elif l == 4:
            colormap.append('orange')  
        elif l == 5:
            colormap.append('purple') 
        elif l == 6: 
            colormap.append('pink')
        elif l == 7: 
            colormap.append('cyan')
    return colormap


# Get Eigenvalues and Eigenvectors of the Laplacian
# Inspired by Daniel Spielman, "Miracles of Algebraic Graph Theory", https://www.youtube.com/watch?v=CDMQR422LGM
def get_eigen(g):
    # First, get the adjacency matrix
    A = nx.adjacency_matrix(g)
    
    # Next generate degrees matrix
    a_shape = A.shape
    a_diagonals = A.sum(axis=1)
    D = scipy.sparse.spdiags(a_diagonals.flatten(),
                             [0],
                             a_shape[0],
                             a_shape[1],
                             format='csr')
    
    # Laplacian
    L = (D - A)
    
    # w are the eigenvalues, sorted in ascending order
    # v are the eigenvectors, v[:,0], v[:,1] etc.
    w, v = eigh(L.todense())
    return w, v


def show_spectralmap(g, nclusters=4):
    w, v = get_eigen(g)
    x = v[:,1]
    y = v[:,2]
    X = [[x[i], y[i]] for i in range(len(x))]
    kmeans = KMeans(n_clusters=nclusters, random_state=0).fit(X)
    
    colormap = label2color(kmeans.labels_)
    ns = list(g.nodes())
    spectral_coordinates = {ns[i] : [x[i], y[i]] for i in range(len(ns))}
    nx.draw(g, pos=spectral_coordinates, node_color=colormap, node_size=[8]*len(ns), with_labels=True) 
    plt.show()


def show_spectralmap_simple(g):
    w, v = get_eigen(g)
    x = v[:,1]
    y = v[:,2]
    ns = list(g.nodes())
    spectral_coordinates = {ns[i] : [x[i], y[i]] for i in range(len(ns))}
    nx.draw(g, pos=spectral_coordinates, node_size=[8]*len(ns), with_labels=True) 
    plt.show()


def show_colormap(g, nclusters=4):
    w, v = get_eigen(g)
    x = v[:,1]
    y = v[:,2]
    X = [[x[i], y[i]] for i in range(len(x))]
    kmeans = KMeans(n_clusters=nclusters, random_state=0).fit(X)
    
    colormap = label2color(kmeans.labels_)
    ns = list(g.nodes())
    nx.draw(g, node_color=colormap, node_size=[8]*len(ns)) 
    plt.show()
