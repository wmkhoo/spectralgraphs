# Spectral graph theory algorithms in Python

Unofficial implementation of [A Spectral Assignment Approach for the Graph Isomorphism Problem](https://arxiv.org/pdf/1411.0969.pdf) by Stefan Klus and Tuhin Sahai

    >>> from spectralgraphs import *

Examples graphs taken from Figures 3c & 3d

    >>> g,h = house_graphs()


Compute vertex mapping based on spectral assignment approach (complexity O(n^3))

    >>> mapping = min_P_la(g,h)
    >>> mapping
    [(1, 2), (2, 1), (5, 6), (7, 7), (6, 5), (3, 4), (4, 3), (8, 8)]


Compute permutation matrix P based on the two-sided orthogonal Procrustes problem (complexity O(2^n))

    >>> P,min = min_P(g,h)
    >>> P
    array([[ 0.40990005, -0.47176473, -0.08882375,  0.08297524,  0.3839172 ,
             0.40990005,  0.52823527, -0.01422941],
           [-0.60799713,  0.40990005, -0.2640839 , -0.08882375,  0.1362324 ,
             0.39200287,  0.40990005,  0.21082674],
           [-0.08882375,  0.08297524,  0.81980009,  0.05647055,  0.3696878 ,
            -0.08882375,  0.08297524,  0.39814661],
           [ 0.1362324 ,  0.3839172 ,  0.34705914,  0.3696878 , -0.36647441,
             0.1362324 ,  0.3839172 , -0.52599811],
           [-0.2640839 , -0.08882375, -0.21599425,  0.81980009,  0.34705914,
            -0.2640839 , -0.08882375, -0.07459434],
           [ 0.40990005,  0.52823527, -0.08882375,  0.08297524,  0.3839172 ,
             0.40990005, -0.47176473, -0.01422941],
           [ 0.39200287,  0.40990005, -0.2640839 , -0.08882375,  0.1362324 ,
            -0.60799713,  0.40990005,  0.21082674],
           [ 0.21082674, -0.01422941, -0.07459434,  0.39814661, -0.52599811,
             0.21082674, -0.01422941,  0.68552182]])
    >>> min
    1.7457077722177437e-14
