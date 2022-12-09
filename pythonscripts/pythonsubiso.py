import networkx as nx
from time import perf_counter
import numpy as np

D = nx.read_edgelist("data/dataEdges.txt", nodetype=int, create_using=nx.DiGraph)
Q = nx.read_edgelist("data/queryEdges.txt", nodetype=int, create_using=nx.DiGraph)
label_dict = {}

fq = open("data/queryLabels.txt", "r")
fq.readline()
for line in fq:
    node, label = line.split()
    node = int(node)
    label = int(label)
    label_dict[node] = label
fq.close()
nx.set_node_attributes(Q, label_dict, "label")

label_dict.clear()
fq = open("data/dataLabels.txt", "r")
fq.readline()
for line in fq:
    node, label = line.split()
    node = int(node)
    label = int(label)
    label_dict[node] = label
fq.close()
nx.set_node_attributes(D, label_dict, "label")

t1_start = perf_counter()
GM = nx.algorithms.isomorphism.DiGraphMatcher(D, Q, node_match=nx.isomorphism.numerical_node_match("label", None))
s = GM.subgraph_monomorphisms_iter()
counter = 0
#fr = open("data/goldenres.txt", "w")
for sub in s:
    #for x in sub:
        #print(x, end=" ", file=fr)
    #print(file=fr)
    #print(sub, file=fr)
    counter = counter + 1
t2_start = perf_counter()
fo = open("data/golden.txt", "w")
print(counter, file=fo)
timeres = t2_start - t1_start
print(str(counter) + " in " + str(timeres))
fo.close()
#fr.close()
