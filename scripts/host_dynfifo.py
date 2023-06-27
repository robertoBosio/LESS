import pynq
from pynq import Overlay
from pynq import allocate
from pprint import pprint
from random import randint as rnd
import networkx as nx
import numpy as np
import time
from time import perf_counter

def test():
    fq = open("/home/xilinx/overlay/data/golden.txt", "r")
    counter_sol_VF2 = int(fq.readline())
    fq.close()

    node_t = np.uint32
    label_t = np.uint8
    datagraph_v = 0
    datagraph_e = 0
    querygraph_v = 0
    querygraph_e = 0

    fq = open("/home/xilinx/overlay/data/datagraph.csv", "r")
    line = fq.readline()
    letter, datagraph_v, datagraph_e = line.split()
    datagraph_v = int(datagraph_v)
    datagraph_e = int(datagraph_e)

    # Allocating space for map (id -> label)
    datagraph_la = np.empty([datagraph_v], dtype=label_t)
    
    # Allocating space for streams. #
    SRC_EDG_D = allocate(shape=(datagraph_e,), dtype=node_t)
    DST_EDG_D = allocate(shape=(datagraph_e,), dtype=node_t)
    SRC_EDG_D_L = allocate(shape=(datagraph_e,), dtype=label_t)
    DST_EDG_D_L = allocate(shape=(datagraph_e,), dtype=label_t)

    for v in range(datagraph_v):
        line = fq.readline()
        letter, node, label, degree = line.split()
        datagraph_la[int(node)] = int(label)
   
    counter = 0
    for e in range(datagraph_e):
        line = fq.readline()
        letter, nodesrc, nodedst = line.split()
        SRC_EDG_D[counter] = int(nodesrc)
        DST_EDG_D[counter] = int(nodedst)
        SRC_EDG_D_L[counter] = int(datagraph_la[int(nodesrc)])
        DST_EDG_D_L[counter] = int(datagraph_la[int(nodedst)])
        counter = counter + 1

    fq.close()

    fq = open("/home/xilinx/overlay/data/querygraph.csv", "r")
    line = fq.readline()
    letter, querygraph_v, querygraph_e = line.split()
    querygraph_v = int(querygraph_v)
    querygraph_e = int(querygraph_e)

    # Allocating space for map (id -> label)
    querygraph_la = np.empty([querygraph_v], dtype=label_t)
    
    # Allocating space for streams. #
    SRC_ORD = allocate(shape=(querygraph_v,), dtype=node_t)
    SRC_EDG_Q = allocate(shape=(querygraph_e,), dtype=node_t)
    DST_EDG_Q = allocate(shape=(querygraph_e,), dtype=node_t)
    SRC_EDG_Q_L = allocate(shape=(querygraph_e,), dtype=label_t)
    DST_EDG_Q_L = allocate(shape=(querygraph_e,), dtype=label_t)

    for v in range(querygraph_v):
        line = fq.readline()
        letter, node, label, degree = line.split()
        querygraph_la[int(node)] = int(label)
   
    counter = 0
    for e in range(querygraph_e):
        line = fq.readline()
        letter, nodesrc, nodedst = line.split()
        SRC_EDG_Q[counter] = int(nodesrc)
        DST_EDG_Q[counter] = int(nodedst)
        SRC_EDG_Q_L[counter] = int(querygraph_la[int(nodesrc)])
        DST_EDG_Q_L[counter] = int(querygraph_la[int(nodedst)])
        counter = counter + 1

    fq.close()

    # Streaming the query order
    fq = open("/home/xilinx/overlay/data/queryOrder.txt", "r")
    counter = 0
    for line in fq:
        SRC_ORD[counter] = int(line)
        counter = counter + 1
    fq.close()
    
    del datagraph_la, querygraph_la
    
    ol = Overlay("/home/xilinx/overlay/design_1.bit")
    mem_counter = 0
    

    FIFO = allocate(shape=(int(40000000/np.dtype(node_t).itemsize),), dtype=node_t)
    MEM = allocate(shape=(int(40000000/np.dtype(node_t).itemsize),), dtype=node_t)

    mem_counter += SRC_ORD.nbytes
    mem_counter += SRC_EDG_Q.nbytes
    mem_counter += DST_EDG_Q.nbytes
    mem_counter += SRC_EDG_Q_L.nbytes
    mem_counter += DST_EDG_Q_L.nbytes
    mem_counter += SRC_EDG_D.nbytes
    mem_counter += DST_EDG_D.nbytes
    mem_counter += SRC_EDG_D_L.nbytes
    mem_counter += DST_EDG_D_L.nbytes
    mem_counter += FIFO.nbytes
    mem_counter += MEM.nbytes
    print("Occupied " + str(mem_counter / (2**20)) + " Mbytes.", flush=True)
   
    SRC_ORD.flush()
    SRC_EDG_Q.flush()
    DST_EDG_Q.flush()
    SRC_EDG_Q_L.flush()
    DST_EDG_Q_L.flush()
    SRC_EDG_D.flush()
    DST_EDG_D.flush()
    SRC_EDG_D_L.flush()
    DST_EDG_D_L.flush()

    ol.subgraphIsomorphism_0.write(0x10, MEM.device_address)
    ol.subgraphIsomorphism_0.write(0x1c, MEM.device_address)
    ol.subgraphIsomorphism_0.write(0x28, MEM.device_address)
    ol.subgraphIsomorphism_0.write(0x34, MEM.device_address)
    ol.subgraphIsomorphism_0.write(0x40, MEM.device_address)
    ol.subgraphIsomorphism_0.write(0x4c, FIFO.device_address)
    
    fq = open("/home/xilinx/overlay/results.txt", "a")
    tot_time_arr = []
    pre_time_arr = []
    #power_arr = []
    #energy_arr = []

    print("Starting", flush=True)
    for x in range(1):
        MEM[:] = 0
        MEM.flush()
        start = perf_counter()

        #Start the kernel
        ol.subgraphIsomorphism_0.write(0x00, 1)

        #First transaction query vertex order
        ol.axi_dma_0.sendchannel.transfer(SRC_ORD)
        ol.axi_dma_0.sendchannel.wait()

        #Second transaction query edges
        ol.axi_dma_3.sendchannel.transfer(DST_EDG_Q_L)
        ol.axi_dma_2.sendchannel.transfer(SRC_EDG_Q_L)
        ol.axi_dma_1.sendchannel.transfer(DST_EDG_Q)
        ol.axi_dma_0.sendchannel.transfer(SRC_EDG_Q)
        ol.axi_dma_0.sendchannel.wait()

        #Third transaction data edges
        ol.axi_dma_3.sendchannel.transfer(DST_EDG_D_L)
        ol.axi_dma_2.sendchannel.transfer(SRC_EDG_D_L)
        ol.axi_dma_1.sendchannel.transfer(DST_EDG_D)
        ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
        ol.axi_dma_0.sendchannel.wait()

        #Fourth transaction data edges
        ol.axi_dma_3.sendchannel.transfer(DST_EDG_D_L)
        ol.axi_dma_2.sendchannel.transfer(SRC_EDG_D_L)
        ol.axi_dma_1.sendchannel.transfer(DST_EDG_D)
        ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
        ol.axi_dma_0.sendchannel.wait()

        #while (not (ol.subgraphIsomorphism_0.read(0x00) & 0x2)):
        #    pass
        while (not(ol.subgraphIsomorphism_0.read(0x58))):
            pass

        end_preprocess = perf_counter()
        #print(ol.subgraphIsomorphism_0.register_map, flush=True)
        while (not(ol.subgraphIsomorphism_0.read(0x70))):
            pass

        end_subiso = perf_counter()
        c = ol.subgraphIsomorphism_0.read(0x68)
        print(ol.subgraphIsomorphism_0.read(0x58))
        
        tot_time = end_subiso - start
        pre_time = end_preprocess - start
        
        #print(ol.subgraphIsomorphism_0.register_map)
        #frame = recorder.frame[recorder.frame['Mark'] == 0]
        #power = frame["5V_power"].mean()
        #energy = power * exec_time
        
        #print(exec_time, exec_time1, exec_time2, exec_time3, sep='\n', file=fq)
        tot_time_arr.append(tot_time)
        pre_time_arr.append(pre_time)
        #power_arr.append(power)
        #energy_arr.append(energy)
        
        if c == counter_sol_VF2:
            print("OK, solutions: ", c, sep="", flush=True)
        else:
            print("NO, expected: ", 
                  counter_sol_VF2, "\tactual: ", c, sep="",  flush=True)
    
    print("tot time: " + str(np.mean(tot_time_arr)) +
            "+-" + str(np.std(tot_time_arr)) + " s", file=fq)
    print("pre time: " + str(np.mean(pre_time_arr)) +
            "+-" + str(np.std(pre_time_arr)) + " s", file=fq)
    #print("power: " + str(np.avg(power_arr)) + "+-" + str(np.std(power_arr)) +
    #        " (energy: " + str(np.avg(energy_arr)) +
    #        "+- " + str(np.std(energy_arr)) + ")")

    print(FIFO[999990:1000000])
    #print(MEM[0:100])
    #print(SRC_EDG_D[:])
    #print(DST_EDG_D[:])
    fq.close()
    del MEM, FIFO, SRC_ORD, SRC_EDG_Q, DST_EDG_Q, SRC_EDG_Q_L, DST_EDG_Q_L, SRC_EDG_D, DST_EDG_D, SRC_EDG_D_L, DST_EDG_D_L
    
for x in range(1):
    test()
