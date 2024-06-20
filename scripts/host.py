import pynq
from pynq import Overlay
from pynq import allocate
from pynq import Clocks
from pynq import PL
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
from time import sleep

def parse_args():

    parser = argparse.ArgumentParser(description='reads ')
    
    parser.add_argument('path', 
            help='input file using absolute path')
    
    args = parser.parse_args()
  
    if not os.path.isdir(args.path):
        print("Wrong overlay directory path.")
        sys.exit(2)
    
    return args

def query_order(querygraph_v, querygraph_e, adjacency_list):
    """ Ordering query node based on degrees. The starting node is 
    the one with highest degree. Then the node with highest number 
    of neighbors in the already ordered set is selected.  """

    # Selecting the node with the highest degree as starting node
    max_degree = 0
    order = []
    query_vertices = list(range(querygraph_v))
    for v in range(querygraph_v):
        degree = len(adjacency_list[v])
        if degree > max_degree:
            max_degree = degree
            start_node = v
    
    order.append(start_node)
    query_vertices.remove(start_node)

    for x in range(querygraph_v - 1):
        max_neigh = 0
        following = query_vertices[0]
        for candidate in query_vertices:
            neighbors_already_matched = 0
            for neighbor in adjacency_list[candidate]:
                if neighbor in order:
                    neighbors_already_matched += 1

            if (neighbors_already_matched > max_neigh):
                max_neigh = neighbors_already_matched
                following = candidate

            # If two nodes have the same number of neighbors already matched, 
            # the one with the highest degree is selected
            elif (neighbors_already_matched == max_neigh):
                if (len(adjacency_list[candidate]) > len(adjacency_list[following])):
                    following = candidate
        
        if max_neigh == 0:
            print("Error: query graph is not connected.")
            return None

        query_vertices.remove(following)
        order.append(following)

    return order

def subiso(test, path):
  
    HASHTABLES_SPACE = 1 << 28  #~ 256 MB
    BLOOM_SPACE  = 1 << 27  #~ 128 MB
    RESULTS_SPACE = 1 << 29  #~ 512 MB
    MAX_QDATA = 300
    BURST_SIZE = 32
    nfile = 0
    tot_time_bench = 0
    time_limit = 2400
    lab_w = 5
    datagraph_v = 0
    datagraph_e = 0
    querygraph_v = 0
    querygraph_e = 0
    mem_counter = 0
    byte_bloom = 16 * 4
    byte_counter = 4
    byte_edge = 5
   
    ## AXILITE register addresses ##
    axi_addresses = {
        "addr_mem0" : 0x10,
        "addr_mem1" : 0x1c,
        "addr_mem2" : 0x28,
        "addr_mem3" : 0x34,
        "addr_bloom" : 0x40,
        "addr_fifo" : 0x4c,
        "addr_qv" : 0x58,
        "addr_qe" : 0x60,
        "addr_de" : 0x68,
        "addr_hash1_w" : 0x74,
        "addr_hash2_w" : 0x7c,
        "addr_dyn_spacel" : 0x84,
        "addr_dyn_spaceh" : 0x88,
        "addr_dyn_ovf" : 0x90,
        "addr_dyn_ovf_ctrl" : 0x94,
        "addr_preproc" : 0xa0,
        "addr_preproc_ctrl" : 0xa4,
        "addr_hit_findminl" : 0x218,
        "addr_hit_findminh" : 0x21c,
        "addr_hit_readmin_counterl" : 0x230,
        "addr_hit_readmin_counterh" : 0x234,
        "addr_hit_readmin_edgel" : 0x248,
        "addr_hit_readmin_edgeh" : 0x24c,
        "addr_hit_intersectl" : 0x260,
        "addr_hit_intersecth" : 0x264,
        "addr_hit_verifyl" : 0x278,
        "addr_hit_verifyh" : 0x27c,
        "addr_req_findminl" : 0x290,
        "addr_req_findminh" : 0x294,
        "addr_req_readmin_counterl" : 0x2a8,
        "addr_req_readmin_counterh" : 0x2ac,
        "addr_req_readmin_edgel" : 0x2c0,
        "addr_req_readmin_edgeh" : 0x2c4,
        "addr_req_intersectl" : 0x2d8,
        "addr_req_intersecth" : 0x2dc,
        "addr_req_verifyl" : 0x2f0,
        "addr_req_verifyh" : 0x2f4,
        "addr_hmsb0_1" : 0x308,
        "addr_hmsb0_2" : 0x30c,
        "addr_hmsb0_ctrl" : 0x310,
        "addr_hmsb1_1" : 0x320,
        "addr_hmsb1_2" : 0x324,
        "addr_hmsb1_ctrl" : 0x328,
        "addr_hmsb2_1" : 0x338,
        "addr_hmsb2_2" : 0x33c,
        "addr_hmsb2_ctrl" : 0x340,
        "addr_hmsb3_1" : 0x350,
        "addr_hmsb3_2" : 0x354,
        "addr_hmsb3_ctrl" : 0x358,
        "addr_resl" : 0x368,
        "addr_resh" : 0x36c,
        "addr_res_ctrl" : 0x370,
    }

    node_t = np.uint32
    label_t = np.uint8
    edge_t = np.dtype([('src', node_t), 
                       ('dst', node_t),
                       ('lsrc', node_t),
                       ('ldst', node_t)]) 

    PL.reset()
    ol = Overlay(path + "design_1.bit", download=False)
    FIFO = allocate(shape=(int(RESULTS_SPACE/np.dtype(edge_t).itemsize),), dtype=edge_t)
    BLOOM = allocate(shape=(BLOOM_SPACE,), dtype=np.uint8)
    MEM = allocate(shape=(int(HASHTABLES_SPACE/np.dtype(node_t).itemsize),), dtype=node_t)
    dynfifo_space = 0

    fres = open(path + "results.txt", "a")
    print("query,datagraph,h1,h2,power,time,preproc,hit_findmin,req_findmin,hit_readmin_counter,req_readmin_counter,hit_readmin_edge,req_readmin_edge,hit_intersect,req_intersect,hit_verify,req_verify,hmsb0,hmsb1,hmsb2,hmsb3,matches", file=fres)

    for data in test.keys():
        datagraph = path + "data/" + data
        counter = 0
        fd = open(datagraph, "r")
        line = fd.readline()
        letter, datagraph_v, datagraph_e = line.split()
        datagraph_v = int(datagraph_v)
        datagraph_e = int(datagraph_e)

        #Computing space for dynamic fifo, aligning it to burst size
        dynfifo_space = datagraph_e + MAX_QDATA;
        dynfifo_space = dynfifo_space - (dynfifo_space % BURST_SIZE) + BURST_SIZE;
        if (dynfifo_space > int(RESULTS_SPACE / np.dtype(edge_t).itemsize)):
            print(f"Error: not enough space for dynamic fifo in {data}.", flush=True)
            continue
        dynfifo_space = int(RESULTS_SPACE / np.dtype(edge_t).itemsize) - dynfifo_space;
        counter = dynfifo_space;

        # Allocating space for map (id -> label)
        datagraph_la = np.empty([datagraph_v], dtype=label_t)

        for v in range(datagraph_v):
            line = fd.readline()
            letter, node, label, degree = line.split()
            datagraph_la[int(node)] = int(label)
        
        print(f"Loading {data} in DDR...", end="", flush=True)
        start = time.perf_counter()

        # Adjust the buffer size for optimal performance
        buffer_size = 4096

        while True:
            lines = fd.readlines(buffer_size)
            if not lines:
                break

            for line in lines:
                letter, nodesrc, nodedst = line.split()
                nodesrc, nodedst = int(nodesrc), int(nodedst)
                FIFO[counter] = (nodesrc, nodedst, datagraph_la[nodesrc], datagraph_la[nodedst])
                counter += 1

        end = time.perf_counter()
        print(f" Done in {end - start:.2f} s", flush=True)
        fd.close()
        del datagraph_la
        
        for querytuple in test[data]:
            ol.download()
            Clocks._instance.PL_CLK_CTRLS[0].DIVISOR0=10

            query = querytuple[0]
            
            querygraph = path + "data/" + query
            counter = dynfifo_space + datagraph_e
            fq = open(querygraph, "r")
            line = fq.readline()
            letter, querygraph_v, querygraph_e = line.split()
            querygraph_v = int(querygraph_v)
            querygraph_e = int(querygraph_e)
    
            # Allocating space for map (id -> label)
            querygraph_la = np.empty([querygraph_v], dtype=label_t)
           
            max_degree = 0
            query_vertices = []
            order = []
            adjacency_list = []
            for v in range(querygraph_v):
                line = fq.readline()
                letter, node, label, degree = line.split()
                querygraph_la[int(node)] = int(label)
                query_vertices.append(int(node))
                degree = int(degree)
                adjacency_list.append([])
                if degree > max_degree:
                    max_degree = degree
                    start_node = int(node)
            
            tablelist = []
            edge_list = []
            for e in range(querygraph_e):
                line = fq.readline()
                letter, nodesrc, nodedst = line.split()
                labelsrc = querygraph_la[int(nodesrc)]
                labeldst = querygraph_la[int(nodedst)]
                nodesrc = int(nodesrc)
                nodedst = int(nodedst)
                adjacency_list[nodesrc].append(nodedst)
                adjacency_list[nodedst].append(nodesrc)
                edge_list.append((nodesrc, nodedst, labelsrc, labeldst))
                
                ## Counting number of tables for memory overflow check
                if (nodesrc < nodedst):
                    direction = True

                tupleedge = (labelsrc, labeldst, direction)

                if tablelist.count(tupleedge) == 0:
                    tablelist.append(tupleedge)

            #Taking as a starting node the one with highest degree 
            order.append(start_node)
            query_vertices.remove(start_node)
            
            for x in range(querygraph_v - 1):
                max_neigh = 0
                for candidate in query_vertices:
                    neighborhood = 0
                    for neighbor in adjacency_list[candidate]:
                        if neighbor in order:
                            neighborhood += 1

                    if (neighborhood > max_neigh):
                        max_neigh = neighborhood
                        following = candidate
                query_vertices.remove(following)
                order.append(following)
            print(f"Order selected before: {order}", flush=True)

            # Reordering the query graph based on the heuristic
            order = query_order(querygraph_v, querygraph_e, adjacency_list)
            #order = [0, 1, 2, 3, 4, 5]
            if order is None:
                continue
            
            print(f"Order selected after: {order}", flush=True)

            # Streaming the query order
            for x in range(querygraph_v):
                FIFO[counter] = (order[x], 0, 0, 0)
                counter = counter + 1
            
            for e in range(querygraph_e):
                FIFO[counter] = edge_list[e]
                counter = counter + 1

            fq.close()
            del querygraph_la

            #Resetting memory space
            MEM.fill(0)
            BLOOM.fill(0)
            MEM.flush()
            BLOOM.flush()
            FIFO.flush()
            
            #hash1_w = int(querytuple[2])
            #hash2_w = int(querytuple[3])

            #Heuristic to compute hash table parameters
            hash1_w = int(0.4 * np.log(5*(10**7)*datagraph_e)) + 2
            hash2_w = int(min(max_degree + 1, 7))

            # Assert that the sum of hash widths is bigger than 14
            if (hash1_w + hash2_w <= 14):
                hash2_w = 14 - hash1_w

            blocks = len(tablelist) * 2**(hash1_w + hash2_w - 14)
            while(blocks > 2048):
                hash2_w -= 1
                blocks = len(tablelist) * 2**(hash1_w + hash2_w - 14)
            
            hashtable_spaceused = len(tablelist) * (2**hash1_w) * (2**hash2_w) * byte_counter
            hashtable_spaceused += datagraph_e * byte_edge
            while (hashtable_spaceused > MEM.nbytes):
                hash2_w -= 1
                hashtable_spaceused = len(tablelist) * (2**hash1_w) * (2**hash2_w) * byte_counter
                hashtable_spaceused += datagraph_e * byte_edge

            bloom_spaceused = len(tablelist) * (2**hash1_w) * byte_bloom
            
            print(f"h1 {hash1_w}, h2 {hash2_w}, blocks: {blocks} max 2048")
            
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_mem0"], MEM.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_mem1"], MEM.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_mem2"], MEM.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_mem3"], MEM.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_bloom"], BLOOM.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_fifo"], FIFO.device_address)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_hash1_w"], hash1_w)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_hash2_w"], hash2_w)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_qv"], querygraph_v)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_qe"], querygraph_e)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_de"], datagraph_e)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_dyn_spacel"], dynfifo_space)
            ol.subgraphIsomorphism_0.write(axi_addresses["addr_dyn_spaceh"], dynfifo_space >> 32)

            mem_counter = 0
            mem_counter += FIFO.nbytes
            mem_counter += MEM.nbytes
            mem_counter += BLOOM.nbytes
            
            #Print useful information on memory occupation
            print(data, querytuple, sep=" ", flush=True)
            print(f"Allocated {(mem_counter / (2**20))} Mb. Hash tables" 
                  f" use {((hashtable_spaceused / MEM.nbytes) * 100):.2f}%,"
                  f" bloom use {((bloom_spaceused / BLOOM.nbytes) * 100):.2f}%."
                  f" {dynfifo_space} lines available in fifo.", flush=True)

            if (hashtable_spaceused <= MEM.nbytes and bloom_spaceused <= BLOOM.nbytes):
                power = []
                start = perf_counter()

                #Start the kernel
                ol.subgraphIsomorphism_0.write(0x00, 1)
                while (not(ol.subgraphIsomorphism_0.read(axi_addresses["addr_preproc_ctrl"]))):
                    pass

                end_preprocess = perf_counter()
                print(f"Preprocessing used {(end_preprocess - start):.3f}s.", flush=True)

# while (not (ol.subgraphIsomorphism_0.read(addr_res_ctrl))):
                while (not (ol.subgraphIsomorphism_0.read(0x00) & 0x2)):
                    with open("/sys/class/hwmon/hwmon2/power1_input") as f_input:
                       power.append(int(f_input.read()))
                    curr_time = perf_counter()
                    if (curr_time - start) > time_limit:
                        print("Failed", flush=True)
                        ol.subgraphIsomorphism_0.write(0x00, 0)
                        break
                    sleep(0.001)
                
                end_subiso = perf_counter()
                overflow = ol.subgraphIsomorphism_0.read(axi_addresses["addr_dyn_ovf"])
                resl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_resl"])
                resh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_resh"])
                if (overflow == 1):
                    print(f"Query not solved due to overflow in partial result fifo.", flush=True)
                else:
                    print(f"Query solved in {(end_subiso - start):.3f}s, using in avg {np.mean(power):.3f}nW. Found {(resh << 32) | resl} matches.", flush=True)
                
                tot_time_bench = tot_time_bench + end_subiso

                hit_findminl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_findminl"])
                hit_findminh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_findminh"])
                hit_readmin_counterl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_readmin_counterl"])
                hit_readmin_counterh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_readmin_counterh"])
                hit_readmin_edgel = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_readmin_edgel"])
                hit_readmin_edgeh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_readmin_edgeh"])
                hit_intersectl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_intersectl"])
                hit_intersecth = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_intersecth"])
                hit_verifyl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_verifyl"])
                hit_verifyh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hit_verifyh"])
                req_findminl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_findminl"])
                req_findminh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_findminh"])
                req_readmin_counterl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_readmin_counterl"])
                req_readmin_counterh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_readmin_counterh"])
                req_readmin_edgel = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_readmin_edgel"])
                req_readmin_edgeh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_readmin_edgeh"])
                req_intersectl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_intersectl"])
                req_intersecth = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_intersecth"])
                req_verifyl = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_verifyl"])
                req_verifyh = ol.subgraphIsomorphism_0.read(axi_addresses["addr_req_verifyh"])
                hmsb0_L = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb0_1"])
                hmsb0_H = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb0_2"])
                hmsb1_L = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb1_1"])
                hmsb1_H = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb1_2"])
                hmsb2_L = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb2_1"])
                hmsb2_H = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb2_2"])
                hmsb3_L = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb3_1"])
                hmsb3_H = ol.subgraphIsomorphism_0.read(axi_addresses["addr_hmsb3_2"])

                if (overflow == 1):
                    resh = 0
                    resl = 0
                    end_subiso = start
                    end_preprocess = start

                print(f"{os.path.basename(querygraph)},",
                      f"{os.path.basename(datagraph)},",
                      f"{hash1_w},{hash2_w}",
                      f",{(np.mean(power)):.3f}",
                      f",{(end_subiso - start):.3f}",
                      f",{(end_preprocess - start):.3f}",
                      f",{(hit_findminh << 32) | hit_findminl}",
                      f",{(req_findminh << 32) | req_findminl}",
                      f",{(hit_readmin_counterh << 32) | hit_readmin_counterl}",
                      f",{(req_readmin_counterh << 32) | req_readmin_counterl}",
                      f",{(hit_readmin_edgeh << 32) | hit_readmin_edgel}",
                      f",{(req_readmin_edgeh << 32) | req_readmin_edgel}",
                      f",{(hit_intersecth << 32) | hit_intersectl}",
                      f",{(req_intersecth << 32) | req_intersectl}",
                      f",{(hit_verifyh << 32) | hit_verifyl}",
                      f",{(req_verifyh << 32) | req_verifyl}",
                      f",{(hmsb0_H << 32) | hmsb0_L}",
                      f",{(hmsb1_H << 32) | hmsb1_L}",
                      f",{(hmsb2_H << 32) | hmsb2_L}",
                      f",{(hmsb3_H << 32) | hmsb3_L}",
                      f",{(resh << 32) | resl}",
                      sep="",
                      file=fres)
                fres.flush()
                #cycles = int((tot_time - pre_time) * 290000000);
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0xb8);
                #empty1 = ol.subgraphIsomorphism_0.read(0xbc);
                #empty = (empty1 << 32) | empty0;
                #print(f"Propose empty         {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0xd0);
                #empty1 = ol.subgraphIsomorphism_0.read(0xb4);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Edgebuild empty       {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0xe8);
                #empty1 = ol.subgraphIsomorphism_0.read(0xec);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 4
                #print(f"Findmin empty         {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x100);
                #empty1 = ol.subgraphIsomorphism_0.read(0x104);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Readmin counter empty {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x118);
                #empty1 = ol.subgraphIsomorphism_0.read(0x11c);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Readmin edge empty    {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x130);
                #empty1 = ol.subgraphIsomorphism_0.read(0x134);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Homomorphism empty    {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x148);
                #empty1 = ol.subgraphIsomorphism_0.read(0x14c);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Batch build empty     {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x160);
                #empty1 = ol.subgraphIsomorphism_0.read(0x164);
                #empty = (empty1 << 32) | empty0;
                #print(f"Tuplebuild empty      {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x178);
                #empty1 = ol.subgraphIsomorphism_0.read(0x17c);
                #empty = (empty1 << 32) | empty0;
                #print(f"Intersect empty       {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x190);
                #empty1 = ol.subgraphIsomorphism_0.read(0x194);
                #empty = (empty1 << 32) | empty0;
                #empty = empty
                #print(f"Bypass filter empty          {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x1a8);
                #empty1 = ol.subgraphIsomorphism_0.read(0x1ac);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Split empty           {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x1c0);
                #empty1 = ol.subgraphIsomorphism_0.read(0x1c4);
                #empty = (empty1 << 32) | empty0;
                #print(f"Verify empty          {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x1d8);
                #empty1 = ol.subgraphIsomorphism_0.read(0x1dc);
                #empty = (empty1 << 32) | empty0;
                #print(f"Compact empty         {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x1f0);
                #empty1 = ol.subgraphIsomorphism_0.read(0x1f4);
                #empty = (empty1 << 32) | empty0;
                #print(f"Filter empty          {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #empty = 0
                #empty0 = ol.subgraphIsomorphism_0.read(0x208);
                #empty1 = ol.subgraphIsomorphism_0.read(0x20c);
                #empty = (empty1 << 32) | empty0;
                #empty = empty * 2
                #print(f"Assembly empty        {empty} cc, {(empty * 100 / cycles):.2f}%", file=fres)
                #print(f"Cycles: {cycles}", file=fres)
#
#                if c == int(querytuple[1]):
#                    print("OK, solutions: ", c, sep="", flush=True)
#                    print(f"{pre_time:.4f}", 
#                          f"{tot_time:.4f}   OK\n",
#                          sep='\n',
#                          file=fres)
#                else:
#                    print("***** NO *****, expected: ", 
#                          querytuple[1], "\tactual: ", c, sep="",  flush=True)
#                    print(f"{pre_time:.4f}", 
#                          f"{tot_time:.4f} **NO**\n",
#                          sep='\n',
#                          file=fres)
                
            else :
                print("Skipped due to memory overflow.", flush=True)
            print(" ", flush=True)

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
    
    # print(f"Total test time: {tot_time_bench:.4f}", file=fres)
    fres.close()
    # del FIFO, GRAPH_SPACE, MEM, BLOOM
    del FIFO, MEM, BLOOM

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
