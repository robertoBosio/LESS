import pynq
from pynq import Overlay
from pynq import allocate
from pprint import pprint
from random import randint as rnd
import networkx as nx
import numpy as np
import time

D = nx.read_edgelist("/home/xilinx/overlay/data/dataEdges.txt", nodetype=int, create_using=nx.DiGraph())
Q = nx.read_edgelist("/home/xilinx/overlay/data/queryEdges.txt", nodetype=int, create_using=nx.DiGraph())
label_dict = {}

fq = open("/home/xilinx/overlay/data/queryLabels.txt", "r")
fq.readline()
for line in fq:
    node, label = line.split()
    node = int(node)
    label = int(label)
    label_dict[node] = label
fq.close()
nx.set_node_attributes(Q, label_dict, "label")

label_dict.clear()
fq = open("/home/xilinx/overlay/data/dataLabels.txt", "r")
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

print("Solutions by VF2: " + str(counter_sol_VF2), flush=True)
numQueryVertex = nx.number_of_nodes(Q)
numDataVertex = nx.number_of_nodes(D)
numQueryEdges = nx.number_of_edges(Q)
numDataEdges = nx.number_of_edges(D)

print("Number of query vertices: " + str(numQueryVertex), flush=True)
print("Number of data vertices: " + str(numDataVertex), flush=True)
print("Number of query edges: " + str(numQueryEdges), flush=True)
print("Number of data edges: " + str(numDataEdges), flush=True)

#16 32-bit word in a 512-bit word
MEM = allocate(shape=(2500*32,), dtype=np.uint16)
MEM1 = allocate(shape=(2500*32,), dtype=np.uint16)
MEM2 = allocate(shape=(2500*32,), dtype=np.uint16)
MEM3 = allocate(shape=(2500*32,), dtype=np.uint16)
MEM4 = allocate(shape=(2500*32,), dtype=np.uint16)
MEM[:] = 0
MEM1[:] = np.load("/home/xilinx/overlay/data/memoryrandom50.npy")
MEM2[:] = np.load("/home/xilinx/overlay/data/memoryrandom50.npy")
MEM3[:] = np.load("/home/xilinx/overlay/data/memoryrandom50.npy")
MEM4[:] = np.load("/home/xilinx/overlay/data/memoryrandom50.npy")

SRC_ORD = allocate(shape=(numQueryVertex,), dtype=np.uint16)

SRC_EDG_Q = allocate(shape=(numQueryEdges,), dtype=np.uint16)
DST_EDG_Q = allocate(shape=(numQueryEdges,), dtype=np.uint16)
SRC_EDG_Q_L = allocate(shape=(numQueryEdges,), dtype=np.uint8)
DST_EDG_Q_L = allocate(shape=(numQueryEdges,), dtype=np.uint8)

SRC_EDG_D = allocate(shape=(numDataEdges,), dtype=np.uint16)
DST_EDG_D = allocate(shape=(numDataEdges,), dtype=np.uint16)
SRC_EDG_D_L = allocate(shape=(numDataEdges,), dtype=np.uint8)
DST_EDG_D_L = allocate(shape=(numDataEdges,), dtype=np.uint8)

RES = allocate(shape=(1000,), dtype=np.uint16)
#RES[:] = 0
ol = Overlay("/home/xilinx/overlay/design_1.bit")
#pprint(ol.ip_dict)
# Streaming the query order
counter = 0
fq = open("/home/xilinx/overlay/data/queryOrder.txt", "r")
for line in fq:
    SRC_ORD[counter] = int(line)
    counter = counter + 1
fq.close()

#print(SRC_ORD)
#print("Written QVO", flush=True)
counter = 0
# Streaming the query edges
for edge in nx.edges(Q):
    SRC_EDG_Q[counter] = edge[0]
    DST_EDG_Q[counter] = edge[1]
    SRC_EDG_Q_L[counter] = Q.nodes[edge[0]]['label']
    DST_EDG_Q_L[counter] = Q.nodes[edge[1]]['label']
    counter = counter + 1


#print(SRC_EDG_Q)
#print(DST_EDG_Q)
#print(SRC_EDG_Q_L)
#print(DST_EDG_Q_L)

#print("Written query edges", flush=True)
counter = 0
# Streaming the query edges
for edge in nx.edges(D):
    SRC_EDG_D[counter] = edge[0]
    DST_EDG_D[counter] = edge[1]
    SRC_EDG_D_L[counter] = D.nodes[edge[0]]['label']
    DST_EDG_D_L[counter] = D.nodes[edge[1]]['label']
    counter = counter + 1

#print(SRC_EDG_D)
#print(DST_EDG_D)
#print(SRC_EDG_D_L)
#print(DST_EDG_D_L)
#print("Written data edges", flush=True)
SRC_ORD.flush()
SRC_EDG_Q.flush()
DST_EDG_Q.flush()
SRC_EDG_Q_L.flush()
DST_EDG_Q_L.flush()
SRC_EDG_D.flush()
DST_EDG_D.flush()
SRC_EDG_D_L.flush()
DST_EDG_D_L.flush()
MEM.flush()
#RES.flush()

#print(MEM[0:100], flush=True)

#print("Addresses axi", flush=True)
ol.subisoWrap_0.write(0x10, MEM.device_address)
ol.subisoWrap_0.write(0x1c, MEM.device_address)
ol.subisoWrap_0.write(0x28, MEM.device_address)
ol.subisoWrap_0.write(0x34, MEM.device_address)
ol.subisoWrap_0.write(0x40, MEM.device_address)

#start = time.time()
#Start the kernel
ol.subisoWrap_0.write(0x00, 1)

#First transaction query vertex order
print("First", flush=True)
ol.axi_dma_0.sendchannel.transfer(SRC_ORD)
ol.axi_dma_0.sendchannel.wait()

#Second transaction query edges
print("Second", flush=True)
ol.axi_dma_2.sendchannel.transfer(DST_EDG_Q_L)
ol.axi_dma_1.sendchannel.transfer(SRC_EDG_Q_L)
ol.axi_dma_3.sendchannel.transfer(DST_EDG_Q)
ol.axi_dma_0.sendchannel.transfer(SRC_EDG_Q)
ol.axi_dma_0.sendchannel.wait()

#Third transaction data edges
print("Third", flush=True)
ol.axi_dma_2.sendchannel.transfer(DST_EDG_D_L)
ol.axi_dma_1.sendchannel.transfer(SRC_EDG_D_L)
ol.axi_dma_3.sendchannel.transfer(DST_EDG_D)
ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
ol.axi_dma_0.sendchannel.wait()

#Fourth transaction data edges
print("Fourth", flush=True)
ol.axi_dma_2.sendchannel.transfer(DST_EDG_D_L)
ol.axi_dma_1.sendchannel.transfer(SRC_EDG_D_L)
ol.axi_dma_3.sendchannel.transfer(DST_EDG_D)
ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
ol.axi_dma_0.sendchannel.wait()

#np.save("/home/xilinx/overlay/memoryrandom50.npy", MEM)
print("0x00: " + str(ol.subisoWrap_0.read(0x00)), flush=True)
ol.axi_dma_4.recvchannel.transfer(RES)

print("dma4 -> 0x32: " + str(ol.axi_dma_4.read(0x34)), flush=True)
print("dma4 -> 0x58: " + str(ol.axi_dma_4.read(0x58)), flush=True)
print("Waiting for result", flush=True)
print(RES.device_address, flush=True)
ol.axi_dma_4.recvchannel.wait()

#end = time.time()
print("End", flush=True)
c = 0
counter = 0
flag = True
while (flag):
    if (RES[c] != 0 and RES[c+1] != 0 and RES[c+2] != 0 and RED[c+3] != 0):
        counter = counter + 1
        c = c + 4
        flag = True
    else:
        flag = False

if counter == counter_sol_VF2:
    print("OK")
else:
    print("NO")

print(MEM[0:100], flush=True)
print(RES[0:10], flush=True)

del MEM, SRC_ORD, SRC_EDG_Q, DST_EDG_Q, SRC_EDG_Q_L, DST_EDG_Q_L, SRC_EDG_D, DST_EDG_D, SRC_EDG_D_L, DST_EDG_D_L, RES
del MEM1, MEM2, MEM3, MEM4







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
