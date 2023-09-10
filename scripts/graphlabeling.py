#!/usr/bin/env python3

import getopt
import argparse
import sys
import os.path 
import numpy as np
import networkx as nx
from random import randint as rnd

def main(): 
    args = parse_args()

    #data_graph = nx.read_edgelist(args.input_file, nodetype=int, create_using=nx.DiGraph)
    #data_graph = nx.read_edgelist(args.input_file, nodetype=int, create_using=nx.Graph)
    vertices_set = set()
    edges = []
    with open(args.input_file, 'r') as input_file:
        for line in input_file:
            parts = line.strip().split()
            if len(parts) == 2:
                vertex1, vertex2 = int(parts[0]), int(parts[1])
                vertices_set.add(vertex1)
                vertices_set.add(vertex2)
                edges.append((vertex1, vertex2))
                
    # Create a set to store all unique vertices
    n_edges = len(edges)
    n_vertices = len(vertices_set)

    labels_edges = np.zeros(n_edges, dtype=np.int32)
    labels_vertices = np.zeros(n_vertices, dtype=np.int32)
    unique_vertices = np.zeros(n_vertices, dtype=np.int32)
    reverse_vertices = np.zeros(n_vertices, dtype=np.int32)
    degree = np.zeros(n_vertices, dtype=np.int32)
    
    remap_vertices_orderedID(unique_vertices, reverse_vertices, degree, args.input_file)

    #Since before 0 is used as flag of not assigned vertex
    for i in range(n_vertices):
        unique_vertices[i] -= 1

    assign_edge_labels(args.labels_edges, n_edges, labels_edges)
    assign_node_labels(args.labels_vertices, n_vertices, labels_vertices)

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
            produce_datagraph_file_DAF(n_vertices, unique_vertices, reverse_vertices, edges, labels_vertices, args.output_path)
        else:
            produce_querygraph_file_DAF(n_vertices, n_edges, unique_vertices, reverse_vertices, edges, labels_vertices, degree, args.output_path)
    if (args.RM):
        print("Producing for RM")
        produce_graph_file_RM(n_vertices, n_edges, unique_vertices, reverse_vertices, edges, labels_vertices, degree, args.output_path)
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
    if (args.GSI):
        print("Producing for GSI")
        produce_graph_file_GSI(n_vertices, n_edges, args.labels_vertices, args.labels_edges, unique_vertices, reverse_vertices, edges, labels_vertices, args.output_path)
    if (args.STMATCH):
        print("Producing for STMatch")
        if (not args.query):
            produce_datagraph_file_STMATCH(edges, unique_vertices, reverse_vertices, labels_vertices, args.output_path)
        else:
            produce_querygraph_file_STMATCH(unique_vertices, reverse_vertices, edges, labels_vertices, args.output_path)
        

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
    parser.add_argument('--GSI',
            action='store_true',
            help='produce for GSI algorithm')
    parser.add_argument('--STMATCH',
            action='store_true',
            help='produce for STMATCH algorithm')
    
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

def remap_vertices_orderedID(unique_vertices, reverse_vertices, degree, input_filename):
    idnum = 1

    #Assgning to each vertex an ID the first time they appear in an edge
    with open(input_filename, 'r') as input_file:
        for line in input_file:
            parts = line.strip().split()
            if len(parts) == 2:
                old_vertex1, old_vertex2 = int(parts[0]), int(parts[1])
                if (unique_vertices[old_vertex1] == 0):
                    unique_vertices[old_vertex1] = idnum
                    reverse_vertices[idnum - 1] = old_vertex1
                    idnum += 1
                if (unique_vertices[old_vertex2] == 0):
                    unique_vertices[old_vertex2] = idnum
                    reverse_vertices[idnum - 1] = old_vertex2
                    idnum += 1
                degree[old_vertex1] += 1
                degree[old_vertex2] += 1

def assign_edge_labels(n_edge_labels, n_edges, edge_labels):
    for i in range(n_edges):
        edge_labels[i] = rnd(0, n_edge_labels)

def assign_node_labels(n_node_labels, n_vertices, vertices_labels):
    for i in range(n_vertices):
        vertices_labels[i] = rnd(0, n_node_labels)

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

def produce_datagraph_file_DAF(n_vertices, vertices, reverse_vertices, edges, labels_vertices, output_path):
    datagraph_file = open(output_path + ".DAF.csv", "w+")
    print("t",
            "1",
            n_vertices, 
            sep = " ",
            file=datagraph_file)
    
    for i, v in enumerate(reverse_vertices):
        print("v", 
        i, 
        labels_vertices[v], 
        sep=" ", 
        file=datagraph_file)

    for edge in edges:
        print("e",
            vertices[edge[0]], 
                vertices[edge[1]],
                "0",
                sep=" ", 
                file=datagraph_file)
    

    datagraph_file.close()

def produce_querygraph_file_DAF(n_vertices, n_edges, vertices, reverse_vertices, edges, labels_vertices, degree, output_path):
    querygraph_file = open(output_path + ".DAF.csv", "w")
    print("t",
            "0",
            n_vertices,
            n_edges * 2,
            sep=" ",
            file=querygraph_file)

    adjacency_lists= []
    for i in range(n_vertices):
        adjacency_lists.append([])

    for edge in edges:
        vertex1 = vertices[edge[0]]
        vertex2 = vertices[edge[1]]
        adjacency_lists[vertex1].append(vertex2)
        adjacency_lists[vertex2].append(vertex1)

    for i, v in enumerate(reverse_vertices):
        neighbors = adjacency_lists[i]
        neighbors.sort()
        print(i, 
                labels_vertices[v],
                degree[v],
                *neighbors,
                sep=" ", 
                file=querygraph_file)
    
    querygraph_file.close()

def produce_graph_file_RM(n_vertices, n_edges, vertices, reverse_vertices, edges, labels_vertices, degree, output_path):
    graph_file = open(output_path + ".RM.csv", "w")
    print("t",
            n_vertices, 
            n_edges, 
            sep = " ",
            file=graph_file)
    
    for i, v in enumerate(reverse_vertices):
        print("v", 
        i, 
        labels_vertices[v], 
        degree[v], 
        sep=" ", 
        file=graph_file)

    for edge in edges:
        print("e",
            vertices[edge[0]], 
                vertices[edge[1]], 
                sep=" ", 
                file=graph_file)
    
    graph_file.close()

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

def produce_datagraph_file_STMATCH(edges, vertices, reverse_vertices, vertices_labels, output_path):
    edges_file = open(output_path + ".STMatch.csv", "w")
    label_file = open(output_path + ".STMatch.label.csv", "w")
    for i, v in enumerate(reverse_vertices):
        print(vertices_labels[v], file=label_file)

    for edge in edges:
        print(vertices[edge[0]], 
                vertices[edge[1]], 
                sep=" ", 
                file=edges_file)
    edges_file.close()
    label_file.close()

def produce_querygraph_file_STMATCH(vertices, reverse_vertices, edges, labels_vertices, output_path):
    graph_file = open(output_path + ".STMatch.csv", "w")
    print("# t 0",
            file=graph_file)
    
    for i, v in enumerate(reverse_vertices):
        print("v", 
        i, 
        labels_vertices[reverse_vertices[i]], 
        sep=" ", 
        file=graph_file)

    for edge in edges:
        print("e",
            vertices[edge[0]], 
                vertices[edge[1]],
                "1", 
                sep=" ", 
                file=graph_file)
    
    graph_file.close()

def produce_graph_file_GSI(n_vertices, n_edges, n_labels_vertices, n_labels_edges, vertices, reverse_vertices, edges, labels_vertices, output_path):
    graph_file = open(output_path + ".GSI.csv", "w")
    print("t # 0",
            file=graph_file)
    
    print(n_vertices, 
            n_edges,
            n_labels_vertices + 1,
            n_labels_edges + 1, 
            sep = " ",
            file=graph_file)
    
    for i, v in enumerate(reverse_vertices):
        print("v", 
        i, 
        labels_vertices[reverse_vertices[i]], 
        sep=" ", 
        file=graph_file)

    for edge in edges:
        print("e",
            vertices[edge[0]], 
                vertices[edge[1]], 
                sep=" ", 
                file=graph_file)
    
    print("t # -1",
            file=graph_file)
    
    graph_file.close()

if __name__ == "__main__":
       main()
