import networkx as nx

if __name__ == '__main__':
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
    counter = 0
    fr = open("data/goldenres.txt", "w")
    for sub in s:
        for x in sub:
            print(x, end=" ", file=fr)
        print(file=fr)
        #print(sub, file=fr)
        counter = counter + 1
    fo = open("data/golden.txt", "w")
    print(counter, file=fo)
    print(counter)
    fo.close()
    fr.close()

