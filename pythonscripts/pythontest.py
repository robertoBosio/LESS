import networkx as nx
import sys, getopt
from time import perf_counter
import numpy as np

def main(argv):
    data_directory = ""
    queries_directory = ""
    n_queries = 1

    try:
        opts, args = getopt.getopt(argv,"d:q:r:")
    except getopt.GetoptError:
        print('pythonlabeling.py -d <data graph directory> -q <queries directory> -r <number of queries>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-q':
            queries_directory = arg
        elif opt == '-d':
            data_directory = arg
        elif opt == '-r':
            n_queries = int(arg)

    for g in range(4):
        fres = open("data/results" + str(g) + ".txt", "w")
        results = np.zeros((n_queries, 10))
        for rep in range(10):
            for q in range(n_queries):
                D = nx.read_edgelist(data_directory + str(g) + "/dataEdges.txt", nodetype=int, create_using=nx.DiGraph)
                Q = nx.read_edgelist(queries_directory + "query" + str(q) +"/queryEdges.txt", nodetype=int, create_using=nx.DiGraph)
                label_dict = {}

                fq = open(queries_directory + "query" + str(q) + "/queryLabels.txt", "r")
                fq.readline()
                for line in fq:
                    node, label = line.split()
                    node = int(node)
                    label = int(label)
                    label_dict[node] = label
                fq.close()
                nx.set_node_attributes(Q, label_dict, "label")

                label_dict.clear()
                fq = open(data_directory + str(g) + "/dataLabels.txt", "r")
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
                #fo = open("data/golden.txt", "w")
                #print(counter, file=fo)
                t2_start = perf_counter()
                timeres = t2_start - t1_start
                print("Query " + str(q) + " in " + str(g) + ": " + str(counter))
                #fo.close()
                #fr.close()
                results[q][rep] = timeres
            print("Graph: " + str(g) + ", rep: " + str(rep))
        print(results, file=fres)
        fres.close()

if __name__ == "__main__":
       main(sys.argv[1:])
