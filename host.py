import pynq
from pynq import Overlay
from pynq import allocate
from random import randint as rnd
import networkx as nx
import numpy as np
import time

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

GM = nx.algorithms.isomorphism.DiGraphMatcher(D, Q, node_match=nx.isomorphism.numerical_node_match("label", None))
s = GM.subgraph_monomorphisms_iter()
counter_sol_VF2 = 0
for sub in s:
    counter_sol_VF2 = counter_sol_VF2 + 1

numQueryVertex = nx.number_of_nodes(Q)
numDataVertex = nx.number_of_nodes(D)
numQueryEdges = nx.number_of_edges(Q)
numDataEdges = nx.number_of_edges(D)

#16 32-bit word in a 512-bit word
MEM = allocate(shape=(2500*16,), dtype=np.uint32)

SRC_ORD = allocate(shape=(numQueryVertex,), dtype=np.uint16)

SRC_EDG_Q = allocate(shape=(numQueryEdges,), dtype=np_uint16)
DST_EDG_Q = allocate(shape=(numQueryEdges,), dtype=np_uint16)
SRC_EDG_Q_L = allocate(shape=(numQueryEdges,), dtype=np_uint8)
DST_EDG_Q_L = allocate(shape=(numQueryEdges,), dtype=np_uint8)

SRC_EDG_D = allocate(shape=(numDataEdges,), dtype=np_uint16)
DST_EDG_D = allocate(shape=(numDataEdges,), dtype=np_uint16)
SRC_EDG_D_L = allocate(shape=(numDataEdges,), dtype=np_uint8)
DST_EDG_D_L = allocate(shape=(numDataEdges,), dtype=np_uint8)

RES = allocate(shape=(100000,), dtype=np.uint16)

ol = Overlay("./subiso_bd/subiso_bd.runs/impl_1/design_1_wrapper.bit")

# Streaming the query order
counter = 0
fq = open("data/queryOrder.txt", "r")
for line in fq:
    SRC_ORD[counter_streams] = int(line))
    counter = counter + 1
fq.close()

counter = 0
# Streaming the query edges
for edge in nx.edges(Q):
    SRC_EDG_Q[counter] = edge[0]
    DST_EDG_Q[counter] = edge[1]
    SRC_EDG_Q_L[counter] = Q.nodes[edge[0]]['label']
    DST_EDG_Q_L[counter] = Q.nodes[edge[1]]['label']
    counter = counter + 1

counter = 0
# Streaming the query edges
for edge in nx.edges(D):
    SRC_EDG_D[counter] = edge[0]
    DST_EDG_D[counter] = edge[1]
    SRC_EDG_D_L[counter] = Q.nodes[edge[0]]['label']
    DST_EDG_D_L[counter] = Q.nodes[edge[1]]['label']
    counter = counter + 1

SRC_ORD.flash()
SRC_EDG_Q.flush()
DST_EDG_Q.flush()
SRC_EDG_Q_L.flush()
DST_EDG_Q_L.flush()
SRC_EDG_D.flush()
DST_EDG_D.flush()
SRC_EDG_D_L.flush()
DST_EDG_D_L.flush()

ol.subisoWrap_0.write(0x10, MEM.device_address)
ol.subisoWrap_0.write(0x1c, MEM.device_address)
ol.subisoWrap_0.write(0x28, MEM.device_address)
ol.subisoWrap_0.write(0x34, MEM.device_address)
ol.subisoWrap_0.write(0x40, MEM.device_address)

start = time.time()
#Start the kernel
ol.subisoWrap_0.write(0x00, 1)

#First transaction query vertex order
ol.axi_dma0.sendchannel.transfer(SRC_ORD)
ol.axi_dma0.sendchannel.wait()

#Second transaction query edges
ol.axi_dma2.sendchannel.transfer(DST_EDG_Q_L)
ol.axi_dma1.sendchannel.transfer(SRC_EDG_Q_L)
ol.axi_dma3.sendchannel.transfer(DST_EDG_Q)
ol.axi_dma0.sendchannel.transfer(SRC_EDG_Q)
ol.axi_dma0.sendchannel.wait()

#Third transaction data edges
ol.axi_dma2.sendchannel.transfer(DST_EDG_D_L)
ol.axi_dma1.sendchannel.transfer(SRC_EDG_D_L)
ol.axi_dma3.sendchannel.transfer(DST_EDG_D)
ol.axi_dma0.sendchannel.transfer(SRC_EDG_D)
ol.axi_dma0.sendchannel.wait()

#Fourth transaction data edges
ol.axi_dma2.sendchannel.transfer(DST_EDG_D_L)
ol.axi_dma1.sendchannel.transfer(SRC_EDG_D_L)
ol.axi_dma3.sendchannel.transfer(DST_EDG_D)
ol.axi_dma0.sendchannel.transfer(SRC_EDG_D)
ol.axi_dma0.sendchannel.wait()

ol.axi_dma0.readchannel.wait()
ol.axi_dma0.readchannel.transfer(RES)

end = time.time()

c = 0
flag = True
while (flag):
    if (RES[c] != 0 and RES[c+1] != 0 and RES[c+2] != and RED[c+3] != 0):
        counter_sol = counter_sol + 1
        c = c + 4
        flag = True
    else:
        flag = False

if counter_sol == counter_sol_VF:
    print("OK")
else:
    print("NO")

del MEM, SRC_ORD, SRC_EDG_Q, DST_EDG_Q, SRC_EDG_Q_L, DST_EDG_Q_L, SRC_EDG_D, DST_EDG_D, SRC_EDG_D_L, DST_EDG_D_L, RES








#data_graph = nx.random_k_out_graph(1000, 4, 3, self_loops=False)
##for edge in nx.edges(data_graph):
##    print(edge[0], edge[1], file=fd) 
##fd.close()
#label_set = ["A", "B", "C", "D", "E", "F", "G", "J", "K", "L", "M"]; 
#label_set[rnd(0, n_labels)]
#
##Assign labels to nodes as attribute
#labels_assigned = {}
#for v in data_graph.nodes():
#    labels_assigned[v] = label_set[rnd(0, 10)]
#nx.set_node_attributes(data_graph, labels_assigned, "label")
