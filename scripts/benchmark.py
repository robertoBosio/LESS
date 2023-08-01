#!/usr/bin/env python3

import getopt
import argparse
import sys
import os.path
import subprocess
import re
from random import randint as rnd
from time import perf_counter
import numpy as np


def main(dir_data, dir_query): 
    #fo = open("/home/roberto/Documents/dataset/benchmark/result.txt", "a")
    fo = open("/home/roberto/Documents/result.txt", "a")
    dir_algo = "/home/roberto/Documents/algorithms"
    
    exec_arrayRM = []
    exec_arrayDAF = []
    output = "i"

    for x in range(10):
        print(x)
        output = subprocess.run([
            dir_algo + "/DAF/bin/daf_parallel_10min",
            "-d", dir_data + "DAF.csv", 
            "-q", dir_query + "DAF.csv",
            "-n", "1",
            "-m", "10000000000",
            "-h", "8"],
            stdout=subprocess.PIPE,
            text=True)
        res = re.search("Average Total Time Per Query : (\d+(\.\d+)) ms", str(output))
        exec_arrayDAF.append(float(res.group(1)))
        
        output = subprocess.run([
            dir_algo + "/RapidMatch/build/matching/RapidMatch.out",
            "-d", dir_data + "RM.csv", 
            "-q", dir_query + "RM.csv",
            "-order", "nd",
            "-time_limit", "300",
            "-num", "MAX"],
            stdout=subprocess.PIPE,
            text=True)
        res = re.search("Query time \(seconds\): (\d+(\.\d+))", str(output))
        exec_arrayRM.append(float(res.group(1)) * 1000)
       
    res = re.search("#Embeddings: ([0-9]+)", str(output))
    print(os.path.basename(dir_data), os.path.basename(dir_query),
            res.group(1), sep=" ", file=fo)
    print("RM: " + str(np.mean(exec_arrayRM)) + "+-" + str(np.std(exec_arrayRM)), file=fo)
    print("DAF: " + str(np.mean(exec_arrayDAF)) + "+-" + str(np.std(exec_arrayDAF)) 
            + "\n", file=fo)
    fo.close()

if __name__ == "__main__":

#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
#
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_facebook",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
#
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/musae_github",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
#    
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/twitter_combined",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
#
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
            "/home/roberto/Documents/dataset/benchmark/queries/query5")
