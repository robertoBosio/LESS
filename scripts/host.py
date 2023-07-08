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
import subprocess
import re
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
# addr_mem3 = 0x40
    addr_bloom = 0x40
    addr_fifo = 0x4c
    addr_qv = 0x58
    addr_qe = 0x60
    addr_de = 0x68
# addr_batch_size = 0x80
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

    ol = Overlay(path + "design_1.bit")

# print(Clocks.fclk0_mhz)
# Clocks._instance.PL_SRC_PLL_CTRLS[0].FBDIV=45
# Clocks._instance.PL_SRC_PLL_CTRLS[0].DIV2=0
# Clocks.fclk0_mhz = 500
# Clocks._instance.PL_SRC_PLL_CTRLS[0].FBDIV=60
    Clocks._instance.PL_CLK_CTRLS[0].DIVISOR0=4
    # Clocks.fclk0_mhz = 250
    #print(Clocks._instance.PL_SRC_PLL_CTRLS)
    #print(Clocks._instance.PL_CLK_CTRLS)
# FIFO = allocate(shape=(int(68000000/np.dtype(node_t).itemsize),), dtype=node_t)
# BLOOM = allocate(shape=((1 << 25),), dtype=np.uint8)
# MEM = allocate(shape=(int((1 << 26)/np.dtype(node_t).itemsize),), dtype=node_t)
    FIFO = allocate(shape=(int(68000000/np.dtype(node_t).itemsize),), dtype=node_t)
    BLOOM = allocate(shape=((1 << 26),), dtype=np.uint8)
    MEM = allocate(shape=(int((1 << 25)/np.dtype(node_t).itemsize),), dtype=node_t)

    fres = open(path + "results.txt", "a")

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
            Clocks._instance.PL_CLK_CTRLS[0].DIVISOR0=4
            print(Clocks.fclk0_mhz)

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
            # order = [ 3, 4, 0, 1, 2]
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
            
            hash1_w = int(querytuple[2])
            hash2_w = int(querytuple[3])
    
            ol.subgraphIsomorphism_0.write(addr_graph, GRAPH_SPACE.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem0, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem1, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_mem2, MEM.device_address)
# ol.subgraphIsomorphism_0.write(addr_mem3, MEM.device_address)
            ol.subgraphIsomorphism_0.write(addr_bloom, BLOOM.device_address)
            ol.subgraphIsomorphism_0.write(addr_fifo, FIFO.device_address)
            ol.subgraphIsomorphism_0.write(addr_hash1_w, hash1_w)
            ol.subgraphIsomorphism_0.write(addr_hash2_w, hash2_w)
            ol.subgraphIsomorphism_0.write(addr_qv, querygraph_v)
            ol.subgraphIsomorphism_0.write(addr_qe, querygraph_e)
            ol.subgraphIsomorphism_0.write(addr_de, datagraph_e)
# ol.subgraphIsomorphism_0.write(addr_batch_size, 2)

           
            hashtable_spaceused = len(tablelist) * (2**hash1_w) * (2**hash2_w) * byte_counter
            hashtable_spaceused += datagraph_e * byte_edge
            bloom_spaceused = len(tablelist) * (2**hash1_w) * byte_bloom

            mem_counter = 0;
            mem_counter += FIFO.nbytes
            mem_counter += MEM.nbytes
            mem_counter += BLOOM.nbytes
            mem_counter += GRAPH_SPACE.nbytes
            
            #Print useful information on memory occupation
            print(data, querytuple, sep=" ", flush=True)
            print(f"Allocated {(mem_counter / (2**20))} Mb. Hash tables" 
                  f" use {((hashtable_spaceused / MEM.nbytes) * 100):.2f}%,"
                  f" bloom use {((bloom_spaceused / BLOOM.nbytes) * 100):.2f}%.", flush=True)

            if (hashtable_spaceused <= MEM.nbytes and bloom_spaceused <= BLOOM.nbytes):
                start = perf_counter()

                #Start the kernel
                ol.subgraphIsomorphism_0.write(0x00, 1)
                while (not(ol.subgraphIsomorphism_0.read(addr_preproc_ctrl))):
                    pass

                end_preprocess = perf_counter()
                print(end_preprocess - start, flush=True)
                checkpoint = end_preprocess

                while (not (ol.subgraphIsomorphism_0.read(0x00) & 0x2)):
                    curr_time = perf_counter()
                    if (curr_time - end_preprocess) > time_limit:
                        print("Failed", flush=True)
                        ol.subgraphIsomorphism_0.write(0x00, 0)
                        break
                    else:
                        if (curr_time - checkpoint) > 10:
                            output = subprocess.run(["xmutil", "platformstats"], 
                                    stdout=subprocess.PIPE,
                                    text=True)
                            res = re.search("([0-9]+) mW", str(output))
                            print(res)
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

                print(f"{os.path.basename(querygraph)}"
                      f" {os.path.basename(datagraph)}"
                      f" h1: {hash1_w} h2: {hash2_w}",
                      file=fres)
                print(ol.subgraphIsomorphism_0.read(addr_dyn_fifo), file=fres)

    # frame = recorder.frame[recorder.frame['Mark'] == 0]
    # power = frame["5V_power"].mean()
    # energy = power * exec_time

                if c == int(querytuple[1]):
                    print("OK, solutions: ", c, sep="", flush=True)
                    print(f"{pre_time:.4f}", 
                          f"{tot_time:.4f}   OK\n",
                          sep='\n',
                          file=fres)
                else:
                    print("***** NO *****, expected: ", 
                          querytuple[1], "\tactual: ", c, sep="",  flush=True)
                    print(f"{pre_time:.4f}", 
                          f"{tot_time:.4f} **NO**\n",
                          sep='\n',
                          file=fres)
                
                ol.download()
            else :
                print("Skipped due to memory overflow.", flush=True)

# nfile += 1
# np.save("/home/ubuntu/" + str(nfile) + ".csv", MEM)
# print("tot time: " + str(np.mean(tot_time_arr)) +
# "+-" + str(np.std(tot_time_arr)) + " s", file=fq)
# print("pre time: " + str(np.mean(pre_time_arr)) +
# "+-" + str(np.std(pre_time_arr)) + " s", file=fq)
# print("power: " + str(np.avg(power_arr)) + "+-" + str(np.std(power_arr)) +
# " (energy: " + str(np.avg(energy_arr)) +
# "+- " + str(np.std(energy_arr)) + ")")
# print(FIFO[4194254:4194354], file=fq)
    
    print(f"Total test time: {tot_time_bench:.4f}", file=fres)
    fres.close()
    del FIFO, GRAPH_SPACE, MEM, BLOOM

if __name__ == "__main__":
    args = parse_args()
    test = {}
    prev_datagraph = ""
    testfile = open(args.path + "test.txt", "r")

    for line in testfile:
        if not(line.startswith("#")):
            datagraph, querygraph, golden, h1, h2 = line.split()
            datagraph = os.path.basename(datagraph)
            querygraph = os.path.basename(querygraph)
            if (datagraph == prev_datagraph):
                test[datagraph].append((querygraph, golden, h1, h2))
            else:
                test[datagraph] = [(querygraph, golden, h1, h2)]
            prev_datagraph = datagraph
    testfile.close()
    
    subiso(test, args.path)
    
    # Allocating space for streams. #
# SRC_EDG_D = allocate(shape=(datagraph_e,), dtype=node_t)
# DST_EDG_D = allocate(shape=(datagraph_e,), dtype=node_t)
# SRC_EDG_D_L = allocate(shape=(datagraph_e,), dtype=label_t)
# DST_EDG_D_L = allocate(shape=(datagraph_e,), dtype=label_t)
    # Allocating space for streams. #
# SRC_ORD = allocate(shape=(querygraph_v,), dtype=node_t)
# SRC_EDG_Q = allocate(shape=(querygraph_e,), dtype=node_t)
# DST_EDG_Q = allocate(shape=(querygraph_e,), dtype=node_t)
# SRC_EDG_Q_L = allocate(shape=(querygraph_e,), dtype=label_t)
# DST_EDG_Q_L = allocate(shape=(querygraph_e,), dtype=label_t)

        #First transaction query vertex order
# ol.axi_dma_0.sendchannel.transfer(SRC_ORD)
# ol.axi_dma_0.sendchannel.wait()

# Second transaction query edges
# ol.axi_dma_3.sendchannel.transfer(DST_EDG_Q_L)
# ol.axi_dma_2.sendchannel.transfer(SRC_EDG_Q_L)
# ol.axi_dma_1.sendchannel.transfer(DST_EDG_Q)
# ol.axi_dma_0.sendchannel.transfer(SRC_EDG_Q)
# ol.axi_dma_0.sendchannel.wait()

# Third transaction data edges
# ol.axi_dma_3.sendchannel.transfer(DST_EDG_D_L)
# ol.axi_dma_2.sendchannel.transfer(SRC_EDG_D_L)
# ol.axi_dma_1.sendchannel.transfer(DST_EDG_D)
# ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
# ol.axi_dma_0.sendchannel.wait()

# Fourth transaction data edges
# ol.axi_dma_3.sendchannel.transfer(DST_EDG_D_L)
# ol.axi_dma_2.sendchannel.transfer(SRC_EDG_D_L)
# ol.axi_dma_1.sendchannel.transfer(DST_EDG_D)
# ol.axi_dma_0.sendchannel.transfer(SRC_EDG_D)
# ol.axi_dma_0.sendchannel.wait()

        #while (not (ol.subgraphIsomorphism_0.read(0x00) & 0x2)):
        #    pass
