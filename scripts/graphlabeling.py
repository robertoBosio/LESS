#!/usr/bin/env python3

import getopt
import argparse
import sys
import os.path 
import networkx as nx
from random import randint as rnd

def main(): 
    args = parse_args()

    #data_graph = nx.read_edgelist(args.input_file, nodetype=int, create_using=nx.DiGraph)
    data_graph = nx.read_edgelist(args.input_file, nodetype=int, create_using=nx.Graph)
    remap_vertices(data_graph)
    assign_edge_labels(data_graph, args.labels_edges)
    assign_node_labels(data_graph, args.labels_vertices)
    if (args.GF):
        print("Producing for GF")
        if (not args.query):
            produce_edge_file_GF(data_graph, args.output_path)
            produce_vertices_file_GF(data_graph, args.output_path)
        else:
            produce_querygraph_file_GF(data_graph, args.output_path)
    if (args.DAF):
        print("Producing for DAF")
        if (not args.query):
            produce_datagraph_file_DAF(data_graph, args.output_path, 0)
        else:
            produce_querygraph_file_DAF(data_graph, args.output_path, 0)
    if (args.RM):
        print("Producing for RM")
        if (not args.query):
            produce_datagraph_file_RM(data_graph, args.output_path)
        else:
            produce_querygraph_file_RM(data_graph, args.output_path)
    if (args.VF3):
        print("Producing for VF3")
        if (not args.query):
            produce_datagraph_file_VF3(data_graph, args.output_path)
        else:
            produce_querygraph_file_VF3(data_graph, args.output_path)
    if (args.CECI):
        print("Producing for CECI")
        if (not args.query):
            produce_datagraph_file_CECI(data_graph, args.output_path)
        else:
            produce_querygraph_file_CECI(data_graph, args.output_path)
        

def parse_args():

    parser = argparse.ArgumentParser(description='reads ')
    
    parser.add_argument('input_file', 
            help='input file using absolute path')
    parser.add_argument('-o',
            '--output_path', 
            default=".",
            help='output directory absolute path.')
    parser.add_argument('-lv', 
            '--labels_vertices', 
            default=2,
            type=int,
            help='number of labels randomly assigned to vertices.')
    parser.add_argument('-le', 
            '--labels_edges', 
            default=1,
            type=int,
            help='number of labels randomly assigned to edges.')
    parser.add_argument('--DAF',
            action='store_true',
            help='produce for DAF algorithm')
    parser.add_argument('--GF',
            action='store_true',
            help='produce for GF algorithm')
    parser.add_argument('--RM',
            action='store_true',
            help='produce for RM algorithm')
    parser.add_argument('--VF3',
            action='store_true',
            help='produce for VF3 algorithm')
    parser.add_argument('--CECI',
            action='store_true',
            help='produce for VF3 algorithm')
    parser.add_argument('--query',
            action='store_true',
            help='produce query')
    
    args = parser.parse_args()
  
    args.labels_vertices = args.labels_vertices - 1
    args.labels_edges = args.labels_edges - 1

    if args.labels_vertices < 0:
        print("Cannot have less than one type of vertex label")
        sys.exit(2)

    if args.labels_vertices < 0:
        print("Cannot have less than one type of edge label")
        sys.exit(2)
    
    if os.path.isdir(args.output_path):
        args.output_path = args.output_path.rstrip("/")
        filename = os.path.splitext(os.path.basename(args.input_file))
        args.output_path = args.output_path + "/" + str(filename[0])
    else:
        print("Wrong output directory path.")
        sys.exit(2)

    if not os.path.isfile(args.input_file):
        print("Wrong input file path.")
        sys.exit(2)

    return args

def remap_vertices(graph):
    node_ids = {}
    n_id = 0
    for node in nx.nodes(graph):
        node_ids[node] = n_id
        n_id = n_id + 1
    nx.set_node_attributes(graph, node_ids, "id")


def assign_edge_labels(graph, n_edge_labels):
    edge_labels = {}
    for edge in nx.edges(graph):
        edge_labels[edge] = rnd(0, n_edge_labels)
    nx.set_edge_attributes(graph, edge_labels, "label")

def assign_node_labels(graph, n_node_labels):
    node_labels = {}
    for node in nx.nodes(graph):
        node_labels[node] = rnd(0, n_node_labels)
    nx.set_node_attributes(graph, node_labels, "label")

def produce_edge_file_GF(graph, output_path):
    edges_file = open(output_path + "GF.csv", "w")
    for edge in nx.edges(graph):
        print(graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                graph[edge[0]][edge[1]]["label"], 
                sep=",", 
                file=edges_file)
    edges_file.close()

def produce_vertices_file_GF(graph, output_path):
    nodes_file = open(output_path + "GF.csv", "w")
    for node in nx.nodes(graph):
        print(graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=",", 
                file=nodes_file)
    nodes_file.close()

def produce_querygraph_file_GF(graph, output_path):
    querygraph_file = open(output_path + "GF.csv", "w")
    listquery = []

    for edge in nx.edges(graph):
        listquery.append("(" + str(graph.nodes[edge[0]]["id"]) + ":" + str(graph.nodes[edge[0]]["label"]) + ")"
                + "-[" + str(graph[edge[0]][edge[1]]["label"]) + "]->"
                + "(" + str(graph.nodes[edge[1]]["id"]) + ":" + str(graph.nodes[edge[1]]["label"]) + ")")

    print(*listquery, sep=", ", file=querygraph_file)
    querygraph_file.close()

def produce_datagraph_file_DAF(graph, output_path, graphid):
    datagraph_file = open(output_path + "DAF.csv", "w+")
    print("t",
            graphid,
            nx.number_of_nodes(graph), 
            sep = " ",
            file=datagraph_file)
    
    for node in nx.nodes(graph):
        print("v", 
                graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=" ", 
                file=datagraph_file)
    
    for edge in nx.edges(graph):
        print("e",
                graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                graph[edge[0]][edge[1]]["label"], 
                sep=" ", 
                file=datagraph_file)

    datagraph_file.close()

def produce_querygraph_file_DAF(graph, output_path, graphid):
    querygraph_file = open(output_path + "DAF.csv", "w")
    print("t",
            graphid,
            nx.number_of_nodes(graph),
            nx.number_of_edges(graph) * 2,
            sep=" ",
            file=querygraph_file)
    
    for node in nx.nodes(graph):
        neighbors = list(nx.neighbors(graph, node))
        neighbors.sort()
        print(graph.nodes[node]["id"], 
                graph.nodes[node]["label"],
                nx.degree(graph, node),
                *neighbors,
                sep=" ", 
                file=querygraph_file)
    
    querygraph_file.close()

def produce_datagraph_file_RM(graph, output_path):
    datagraph_file = open(output_path + "RM.csv", "w")
    print("t",
            nx.number_of_nodes(graph), 
            nx.number_of_edges(graph), 
            sep = " ",
            file=datagraph_file)
    
    for node in nx.nodes(graph):
        print("v", 
                graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                nx.degree(graph, node),
                sep=" ", 
                file=datagraph_file)
    
    for edge in nx.edges(graph):
        print("e",
                graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                #graph[edge[0]][edge[1]]["label"], 
                sep=" ", 
                file=datagraph_file)

    datagraph_file.close()

def produce_querygraph_file_RM(graph, output_path):
    querygraph_file = open(output_path + "RM.csv", "w")
    print("t",
            nx.number_of_nodes(graph), 
            nx.number_of_edges(graph), 
            sep = " ",
            file=querygraph_file)
    
    for node in nx.nodes(graph):
        print("v", 
                graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                nx.degree(graph, node),
                sep=" ", 
                file=querygraph_file)
    
    for edge in nx.edges(graph):
        print("e",
                graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                #graph[edge[0]][edge[1]]["label"], 
                sep=" ", 
                file=querygraph_file)

    querygraph_file.close()

def produce_datagraph_file_VF3(graph, output_path):
    datagraph_file = open(output_path + "VF3.csv", "w")
    print("#Number of nodes",
            file=datagraph_file)

    print(nx.number_of_nodes(graph), 
            file=datagraph_file)
    
    print("\n#Node attributes",
            file=datagraph_file)

    for node in nx.nodes(graph):
        print(graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=" ", 
                file=datagraph_file)
    
    for node in nx.nodes(graph):
        print("\n#edges of node",
                node,
                sep=" ",
                file=datagraph_file)

        print(nx.degree(graph, node),
                file=datagraph_file)

        for neighbor in nx.neighbors(graph, node):
            print(graph.nodes[node]["id"], 
                    graph.nodes[neighbor]["id"], 
                    graph[node][neighbor]["label"], 
                    sep=" ", 
                    file=datagraph_file)

    datagraph_file.close()

def produce_querygraph_file_VF3(graph, output_path):
    querygraph_file = open(output_path + "VF3.csv", "w")
    print("#Number of nodes",
            file=querygraph_file)

    print(nx.number_of_nodes(graph), 
            file=querygraph_file)
    
    print("\n#Node attributes",
            file=querygraph_file)
    
    for node in nx.nodes(graph):
        print(graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=" ", 
                file=querygraph_file)

    for node in nx.nodes(graph):
        print("\n#edges of node",
                node,
                sep=" ",
                file=querygraph_file)

        print(nx.degree(graph, node),
                file=querygraph_file)

        for neighbor in nx.neighbors(graph, node):
            print(graph.nodes[node]["id"], 
                    graph.nodes[neighbor]["id"], 
                    graph[node][neighbor]["label"], 
                    sep=" ", 
                    file=querygraph_file)

    querygraph_file.close()

def produce_datagraph_file_CECI(graph, output_path):
    datagraph_file = open(output_path + "CECI.csv", "w")
    print("t # 0",
            file=datagraph_file)
    
    for node in nx.nodes(graph):
        print("v", 
                graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=" ", 
                file=datagraph_file)
    
    for edge in nx.edges(graph):
        print("e",
                graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                graph[edge[0]][edge[1]]["label"], 
                sep=" ", 
                file=datagraph_file)

    datagraph_file.close()

def produce_querygraph_file_CECI(graph, output_path):
    querygraph_file = open(output_path + "CECI.csv", "w")
    print("t # 0",
            file=querygraph_file)
    
    for node in nx.nodes(graph):
        print("v", 
                graph.nodes[node]["id"], 
                graph.nodes[node]["label"], 
                sep=" ", 
                file=querygraph_file)
    
    for edge in nx.edges(graph):
        print("e",
                graph.nodes[edge[0]]["id"], 
                graph.nodes[edge[1]]["id"], 
                graph[edge[0]][edge[1]]["label"], 
                sep=" ", 
                file=querygraph_file)

    querygraph_file.close()

if __name__ == "__main__":
       main()
