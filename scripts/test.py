#!/usr/bin/env python3

import getopt
import argparse
import sys
import os.path
import subprocess
import re
from random import randint as rnd

def main(): 
    args = parse_args()
    memory_avail = 2**27
    dir_algo = "/home/user/Documents/algorithms"
    dir_data = (str(args.path_data) + str(args.input_file_data))
    dir_query = (str(args.path_query) + str(args.input_file_query))

    if not os.path.isfile(dir_data + "RM.csv"):
        print("Wrong data graph file path.")
        sys.exit(2)

    if not os.path.isfile(dir_query + "RM.csv"):
        print("Wrong query graph file path.")
        sys.exit(2)


    datagraph_file = open(dir_data + "RM.csv", "r")
    line = datagraph_file.readline().split()
    datagraph_vertices = int(line[1])
    datagraph_edges = int(line[2])
    datagraph_file.close()

    querygraph_file = open(dir_query + "RM.csv", "r")
    line = querygraph_file.readline().split()
    querygraph_vertices = int(line[1])
    querygraph_edges = int(line[2])
    querygraph_file.close()

    output = subprocess.run([
        dir_algo + "/RapidMatch/build/matching/RapidMatch.out",
        "-d", dir_data + "RM.csv", 
        "-q", dir_query + "RM.csv",
        "-order", "nd",
        "-num", "MAX"],
            stdout=subprocess.PIPE,
            text=True)
       
    res = re.search("#Embeddings: ([0-9]+)", str(output))

    print("Data graph -> vertices: ", datagraph_vertices,
            " edges: ", datagraph_edges, 
            "\nQuery graph -> vertices: ", querygraph_vertices, 
            " edges: ", querygraph_edges,
            sep = "")
    print("\n", res.group(0), sep="")
    
    kernel_used =  ((2**(args.hash1 + args.hash2) * args.counter_width) * querygraph_edges
            + datagraph_edges * args.vertex_size * 2)

    res_used = int(res.group(1)) * querygraph_vertices * args.vertex_size

    dma_used = (querygraph_vertices * args.vertex_size + 
            querygraph_edges * (2 * args.vertex_size + 2) +
            datagraph_edges * (2 * args.vertex_size + 2))

    bytes_used = kernel_used + res_used + dma_used

    print("Used ", int(bytes_used), " bytes, ", 
            (int(bytes_used)*100/128000000), "% of 128Mb\n",
            "\t", kernel_used, " bytes for tables, ", 
            (int(kernel_used)*100/128000000), "% of 128Mb\n",
            "\t", dma_used, " bytes for dma, ", 
            (int(dma_used)*100/128000000), "% of 128Mb\n",
            "\t", res_used, " bytes for results, ", 
            (int(res_used)*100/128000000), "% of 128Mb\n",
            sep="")
    
    if ((memory_avail - bytes_used) < res_used):
        print("-- Space might not be enough --")


def parse_args():

    parser = argparse.ArgumentParser(description='reads ')
    
    parser.add_argument('input_file_data', 
            help='data graph input file using absolute path')
    parser.add_argument('input_file_query', 
            help='query graph input file using absolute path')
    parser.add_argument('-pq',
            '--path_query',
            help='directory path of the query',
            default="/home/user/Documents/dataset/benchmark/queries/")
    parser.add_argument('-pd',
            '--path_data',
            help='directory path of the data',
            default="/home/user/Documents/dataset/benchmark/labelled/")
    parser.add_argument('-s',
            '--vertex_size',
            help='size of vertex id in bytes',
            default=2,
            type=int)
    parser.add_argument('-h1',
            '--hash1',
            help='size of hash1 in bits',
            default=5,
            type=int)
    parser.add_argument('-h2',
            '--hash2',
            help='size of hash2 in bits',
            default=5,
            type=int)
    parser.add_argument('-c',
            '--counter_width',
            help='size of counter in bytes',
            default=2,
            type=int)
    
    args = parser.parse_args()
  
    return args

if __name__ == "__main__":
       main()
