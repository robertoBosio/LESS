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
    fo = open("/home/roberto/Documents/result_benchmark.csv", "a")
    dir_algo = "/home/roberto/Documents/algorithms"
    
    exec_arrayRM = []
    watt_arrayRM = []
    exec_arrayDAF = []
    watt_arrayDAF = []
    
    output = "i"
    print(f"{dir_data} {dir_query}")
    for x in range(5):
        print(x)
        command = f"sudo perf stat -e power/energy-ram/,power/energy-pkg/ -v {dir_algo}/DAF/bin/daf_parallel_10min -d {dir_data}DAF.csv -q {dir_query}DAF.csv -n 1 -m 1000000000 -h 1" 
        output = subprocess.run(command, shell=True, capture_output=True, text=True)
        match_time = re.search("Average Total Time Per Query : (\d+(\.\d+)) ms", str(output.stdout))
        resDAF = re.search("Total Number of Found Matches: (\d{1,3}(?:,\d{3})*)", str(output.stdout))
        power_ram = re.search("(\d+(\.\d+)) Joules power/energy-ram/", str(output.stderr))
        power_cpu = re.search("(\d+(\.\d+)) Joules power/energy-pkg/", str(output.stderr))
        tot_time = re.search("(\d+(\.\d+)) seconds time elapsed", str(output.stderr))
        resDAF = str(resDAF.group(1))
        resDAF = resDAF.replace(',','')
        resDAF = int(resDAF)
        daf_time = float(match_time.group(1))
        joule_ram = float(power_ram.group(1))
        joule_cpu = float(power_cpu.group(1))
        perf_time = float(tot_time.group(1))
        watt_tot = ((joule_cpu + joule_ram) / perf_time)
        daf_time = daf_time / 1000
        exec_arrayDAF.append(daf_time)
        watt_arrayDAF.append(watt_tot)
        
        command = f"sudo perf stat -e power/energy-ram/,power/energy-pkg/ -v {dir_algo}/RapidMatch/build/matching/RapidMatch.out -d {dir_data}RM.csv -q {dir_query}RM.csv -order nd -time_limit 300 -num MAX" 
        output = subprocess.run(command, shell=True, capture_output=True, text=True)
        match_time = re.search("Query time \(seconds\): (\d+(\.\d+))", str(output.stdout))
        power_ram = re.search("(\d+(\.\d+)) Joules power/energy-ram/", str(output.stderr))
        power_cpu = re.search("(\d+(\.\d+)) Joules power/energy-pkg/", str(output.stderr))
        tot_time = re.search("(\d+(\.\d+)) seconds time elapsed", str(output.stderr))
        resRM = re.search("#Embeddings: ([0-9]+)", str(output.stdout))
        resRM = int(resRM.group(1))
        rm_time = float(match_time.group(1))
        joule_ram = float(power_ram.group(1))
        joule_cpu = float(power_cpu.group(1))
        perf_time = float(tot_time.group(1))
        watt_tot = ((joule_cpu + joule_ram) / perf_time)
        exec_arrayRM.append(rm_time)
        watt_arrayRM.append(watt_tot)
       
#    print(os.path.basename(dir_data), os.path.basename(dir_query),
#            res.group(1), sep=" ", file=fo)
#    print("RM: " + str(np.mean(exec_arrayRM)) + " " + str(np.std(exec_arrayRM)), file=fo)
#    print("DAF: " + str(np.mean(exec_arrayDAF)) + " " + str(np.std(exec_arrayDAF)) 
#            + "\n", file=fo)
    print(f"RM\t{resRM}\t{(np.mean(exec_arrayRM)):.3f}\t{(np.std(exec_arrayRM)):.3f}\t{(np.mean(watt_arrayRM)):.3f}\t{(np.std(watt_arrayRM)):.3f}", file=fo)
    print(f"DAF\t{resDAF}\t{(np.mean(exec_arrayDAF)):.3f}\t{(np.std(exec_arrayDAF)):.3f}\t{(np.mean(watt_arrayDAF)):.3f}\t{(np.std(watt_arrayDAF)):.3f}", file=fo)
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
    main("/home/roberto/Documents/dataset/benchmark/labelled/email-Enron",
            "/home/roberto/Documents/dataset/benchmark/queries/query4")
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
#    main("/home/roberto/Documents/dataset/benchmark/labelled/dblp",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
#
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query0")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query1")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query2")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query3")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query4")
#    main("/home/roberto/Documents/dataset/benchmark/labelled/youtube",
#            "/home/roberto/Documents/dataset/benchmark/queries/query5")
