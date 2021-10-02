from spectralgraphs import *

# Graphs taken from Figures 3c & 3d, Stefan Klus and Tuhin Sahai, "A Spectral Assignment Approach for the Graph Isomorphism Problem", https://arxiv.org/pdf/1411.0969.pdf
g, h = house_graphs()

P, min = min_P(g,h)
print(P)
print(min)

mapping = min_P_la(g,h)
print(mapping)