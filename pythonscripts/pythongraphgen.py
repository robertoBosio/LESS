from random import randint as rnd
import networkx as nx

if __name__ == '__main__':
    fd = open("data/dataEdges.txt", "w")
    data_graph = nx.random_k_out_graph(1000, 4, 3, self_loops=False)
    print(nx.number_of_edges(data_graph), file=fd)
    for edge in nx.edges(data_graph):
        print(edge[0], edge[1], file=fd) 
    fd.close()
