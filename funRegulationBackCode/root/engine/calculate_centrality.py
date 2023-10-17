import networkx as nx
from networkx.algorithms.centrality import degree_centrality,closeness_centrality, betweenness_centrality, eigenvector_centrality, harmonic_centrality


### Cria o grafo baseado na lib neworkX usando as
### infos recebidas com os genes ltf e tg do fungo
def create_graph(tf_nodes, tg_nodes):                
    edges = list(zip(tf_nodes, tg_nodes))
    DG = nx.DiGraph()
    DG.add_edges_from(edges)
    return DG

def calculate_degree_centrality(graph):
    return degree_centrality(graph)

def calculate_closeness_centrality(graph):
    return closeness_centrality(graph)

def calculate_betweenness_centrality(graph):
    return betweenness_centrality(graph)

def calculate_eigenvector_centrality(graph):
    return eigenvector_centrality(graph)

def calculate_harmonic_centrality(graph):
    return harmonic_centrality(graph)

