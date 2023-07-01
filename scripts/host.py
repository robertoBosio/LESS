import pynq
from pynq import Overlay
from pynq import allocate
from pynq import Clocks
from pprint import pprint
from random import randint as rnd
import numpy as np
import time
import getopt
import argparse
import os.path
from time import perf_counter

def parse_args():

    parser = argparse.ArgumentParser(description='reads ')
    
    parser.add_argument('path', 
            help='input file using absolute path')
    
    args = parser.parse_args()
  
    if not os.path.isdir(args.path):
        print("Wrong overlay directory path.")
        sys.exit(2)
    
    return args

def subiso(test, path):
  
    nfile = 0
    tot_time_bench = 0
    time_limit = 4000
    lab_w = 5
    datagraph_v = 0
    datagraph_e = 0
    querygraph_v = 0
    querygraph_e = 0
    mem_counter = 0
    byte_bloom = 16 * 8
    byte_counter = 4
    byte_edge = 5
   
    ## AXILITE register addresses ##
    addr_graph = 0x10
    addr_mem0 = 0x1c
    addr_mem1 = 0x28
    addr_mem2 = 0x34
    addr_bloom = 0x40
    addr_fifo = 0x4c
    addr_qv = 0x58
    addr_qe = 0x60
    addr_de = 0x68
    addr_hash1_w = 0x74
    addr_hash2_w = 0x7c
    addr_dyn_fifo = 0x84
    addr_dyn_fifo_ctrl = 0x8c
    addr_preproc = 0x9c
    addr_preproc_ctrl = 0xa0
    addr_res = 0xac
    addr_res_ctrl = 0xb4

    node_t = np.uint32
    label_t = np.uint8
    edge_t = np.dtype([('src', node_t), 
                       ('dst', node_t)]) 

    # Load the overlay
    ol = Overlay(path + "design_1.bit")

    # Allocate memory buffers
    FIFO = allocate(shape=(int(68000000/np.dtype(node_t).itemsize),), dtype=node_t)
    BLOOM = allocate(shape=((1 << 26),), dtype=np.uint8)
    MEM = allocate(shape=(int((1 << 25)/np.dtype(node_t).itemsize),), dtype=node_t)

    for data in test.keys():
        datagraph = path + "data/" + data
        counter = 0
        fd = open(datagraph, "r")
        line = fd.readline()
        letter, datagraph_v, datagraph_e = line.split()
        datagraph_v = int(datagraph_v)
        datagraph_e = int(datagraph_e)

        # Allocating space for map (id -> label)
        datagraph_la = np.empty([datagraph_v], dtype=label_t)

        for v in range(datagraph_v):
            line = fd.readline()
            letter, node, label, degree = line.split()
            datagraph_la[int(node)] = int(label)
    
        GRAPH_SPACE = allocate(shape=(datagraph_e + 100,), dtype=edge_t)

        print("Loading", data, "in DDR...", sep=" ", end = "", flush=True)
        start = perf_counter()
        for e in range(datagraph_e):
            line = fd.readline()
            letter, nodesrc, nodedst = line.split()
            nodesrc = (int(nodesrc) << lab_w) | datagraph_la[int(nodesrc)]
            nodedst = (int(nodedst) << lab_w) | datagraph_la[int(nodedst)]
            GRAPH_SPACE[counter] = (nodesrc, nodedst)
            counter = counter + 1
        end = perf_counter()
        print(" Done in ", end - start, "s", sep="", flush=True)

        fd.close()
        del datagraph_la
        
        for querytuple in test[data]:

            # Force the clock at 250MHz given the PLL at 1GHz
            # There is a discrepancy between the PLL set by Vivado
            # and the real one on the KRIA
            Clocks._instance.PL_CLK_CTRLS[0].DIVISOR0=4

            query = querytuple[0]
            
            querygraph = path + "data/" + query
            counter = datagraph_e
            fq = open(querygraph, "r")
            line = fq.readline()
            letter, querygraph_v, querygraph_e = line.split()
            querygraph_v = int(querygraph_v)
            querygraph_e = int(querygraph_e)
    
            # Allocating space for map (id -> label)
            querygraph_la = np.empty([querygraph_v], dtype=label_t)
           
            for v in range(querygraph_v):
                line = fq.readline()
                letter, node, label, degree = line.split()
                querygraph_la[int(node)] = int(label)
            
            # Streaming the query order
            for x in range(querygraph_v):
                GRAPH_SPACE[counter] = (int(x), 0)
                counter = counter + 1
            
            tablelist = []
            for e in range(querygraph_e):
                line = fq.readline()
                letter, nodesrc, nodedst = line.split()
                labelsrc = querygraph_la[int(nodesrc)] 
                labeldst = querygraph_la[int(nodedst)] 
                nodesrc = (int(nodesrc) << lab_w) | labelsrc
                nodedst = (int(nodedst) << lab_w) | labeldst
                GRAPH_SPACE[counter] = (nodesrc, nodedst)
                counter = counter + 1
                
                ## Counting number of tables for memory overflow check
                if (nodesrc < nodedst):
                    direction = True

                tupleedge = (labelsrc, labeldst, direction)

                if tablelist.count(tupleedge) == 0:
                    tablelist.append(tupleedge)

            fq.close()
            del querygraph_la

            #Resetting memory space
            MEM.fill(0)
            BLOOM.fill(0)
            MEM.flush()
            BLOOM.flush()
            GRAPH_SPACE.flush()
            
            hash1_w = 12
            hash2_w = 5
    
            ol.subgraphIsomorphism_0.write(addr_graph, GRAPH_SPACE.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem0, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem1, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem2, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_bloom, BLOOM.device_address)
            ol.subgraphIsomorphism_0.write(addr_fifo, FIFO.device_address)
            ol.subgraphIsomorphism_0.write(addr_hash1_w, hash1_w)
            ol.subgraphIsomorphism_0.write(addr_hash2_w, hash2_w)
            ol.subgraphIsomorphism_0.write(addr_qv, querygraph_v)
            ol.subgraphIsomorphism_0.write(addr_qe, querygraph_e)
            ol.subgraphIsomorphism_0.write(addr_de, datagraph_e)
           
            hashtable_spaceused = len(tablelist) * (2**hash1_w) * (2**hash2_w) * byte_counter
            hashtable_spaceused += datagraph_e * byte_edge
            bloom_spaceused = len(tablelist) * (2**hash1_w) * byte_bloom

            mem_counter = 0;
            mem_counter += FIFO.nbytes
            mem_counter += MEM.nbytes
            mem_counter += BLOOM.nbytes
            mem_counter += GRAPH_SPACE.nbytes

            if (hashtable_spaceused <= MEM.nbytes and bloom_spaceused <= BLOOM.nbytes):
                start = perf_counter()

                #Start the kernel
                ol.subgraphIsomorphism_0.write(0x00, 1)
                while (not(ol.subgraphIsomorphism_0.read(addr_preproc_ctrl))):
                    pass

                end_preprocess = perf_counter()
                checkpoint = end_preprocess

                while (not (ol.subgraphIsomorphism_0.read(0x00) & 0x2)):
                    curr_time = perf_counter()
                    if (curr_time - end_preprocess) > time_limit:
                        print("Failed.", flush=True)
                        ol.subgraphIsomorphism_0.write(0x00, 0)
                        break
                    else:
                        if (curr_time - checkpoint) > 10:
                            print(ol.subgraphIsomorphism_0.read(addr_dyn_fifo), ", ", curr_time - end_preprocess, "s", sep="", flush=True)
                            checkpoint = curr_time
                        else :
                            pass

                end_subiso = perf_counter()
                c = ol.subgraphIsomorphism_0.read(addr_res)

                tot_time = end_subiso - start
                pre_time = end_preprocess - start
                tot_time_bench = tot_time_bench + tot_time

                if (tot_time > time_limit):
                    tot_time = "failed"

                if c == int(querytuple[1]):
                    print(f"OK: expected {querytuple[1]}, \tfound: {c}", flush=True)
                else:
                    print(f"***** NO *****, expected {querytuple[1]}, \tfound {c}",flush=True)
                
                ol.download()
            else :
                print("Skipped due to memory overflow.", flush=True)

    del FIFO, GRAPH_SPACE, MEM, BLOOM

if __name__ == "__main__":
    args = parse_args()
    test = {}
    prev_datagraph = ""
    testfile = open(args.path + "run_list.txt", "r")

    for line in testfile:
        if not(line.startswith("#")):
            datagraph, querygraph, golden = line.split()
            datagraph = os.path.basename(datagraph)
            querygraph = os.path.basename(querygraph)
            if (datagraph == prev_datagraph):
                test[datagraph].append((querygraph, golden))
            else:
                test[datagraph] = [(querygraph, golden)]
            prev_datagraph = datagraph
    testfile.close()
    
    subiso(test, args.path)
    
