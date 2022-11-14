from random import randint as rnd
import getopt, sys

def main(argv):
    n_labels = 2
    try:
        opts, args = getopt.getopt(argv,"l:")
    except getopt.GetoptError:
        print('pythonlabeling.py -l <number of labels>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-l':
            if int(arg) > 0 and int(arg) < 12:
                n_labels = int(arg)
    
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
        print(node, ord(label_set[rnd(0, n_labels)]), sep=" ", file=fl)

    fd.close()

if __name__ == "__main__":
       main(sys.argv[1:])
