from random import randint as rnd
if __name__ == '__main__':
    fd = open("data/dataEdges.txt", "r")
    label_set = ["A", "B", "C", "D", "E", "F", "G", "J", "K", "L", "M"]; 
    set_nodes = set()
    fd.readline()
    for line in fd:
        nodesrc, nodedst = line.split()
        nodesrc = int(nodesrc)
        nodedst = int(nodedst)
        set_nodes.add(nodesrc)
        set_nodes.add(nodedst)

    fl = open("data/dataLabels.txt", "w")
    print(len(set_nodes), file=fl)
    for node in set_nodes:
        print(node, ord(label_set[rnd(0, 5)]), sep=" ", file=fl)

    fd.close()
