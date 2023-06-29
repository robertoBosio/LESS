# Subgraph isomorphism acceleration on FPGAs using High-Level Synthesis 

[![vitis_hls](https://img.shields.io/badge/vitis--hls-2022.2-blue)](https://docs.xilinx.com/r/2022.2-English/ug1399-vitis-hls/Introduction)
[![vivado](https://img.shields.io/badge/vivado-2022.2-blue)](https://docs.xilinx.com/r/en-US/ug901-vivado-synthesis/Vivado-Documentation)
[![license](https://img.shields.io/badge/license-BSD--3--Clause%20-blue)](https://github.com/robertoBosio/subgraph-iso/blob/master/LICENSE)

## Introduction

Subgraph isomorphism (or subgraph matching) is a well-known NP-hard problem that consists of searching all the distinct embeddings of a query graph in a large data graph. It has a wide range of applications, almost in all the domains in which graph patterns reveal valuable information, especially in social network analysis, chemical compound search, and computer-aided design.
The project aims to shows a all-in-FPGA solution devolped using an High-Level Synthesis approach.

## Getting Started

To compile and upload the kernel on a Xilinx KRIA™ KV26, please follow the instructions below. This section will guide you through the necessary steps to set up the project and deploy it on the FPGA.

### Prerequisites

Before running the project, ensure that you have the following prerequisites installed:

- Vivado 2022.2, with KRIA™ KV26 board files installed.
- Vitis HLS 2022.2
- Pynq and NumPy installed on KRIA™ KV26.

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/robertoBosio/subgraph-iso.git
   ```

2. Navigate to the project directory:

   ```bash
   cd subgraph-iso
   ```

3. Setup the paths:

   ```bash
   export XILINX_VIVADO=/path/to/Vivado/directory
   export XILINX_HLS=/path/to/Vitis_HLS/directory
   ```
4. Compile, implement and generate bitstream:

   ```bash
   make run IMPL=1
   ```

### Input file format
Both the input query graph and data graph are vertex-labeled. Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. The VertexId should start from 0 and be in the range [0, N - 1]. The following is a correct input sample. Other examples are stored under the dataset folder.

Example:
```bash
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

### Run

1. To search a query inside a graph add an entry formatted as '/path/to/datagraph /path/to/query #solutions' in the file run_list.txt which is located under the scripts  directory.

    Example:
    ```bash
    ../dataset/graphs/email-Enron.csv ../dataset/queries/query0.csv 38152
    ```

2. Change the KRIA's IP address in the file script_load.sh under the scripts directory.

    ```bash
        user=root
        path=/home/ubuntu
    ->  ip=192.168.99.170
    ```
3. Make the script executable by running:
    ```bash
    chmod +x script_load.sh
    ```
4. Launch the script by running:
    ```bash
    ./script_load.sh
    ```
The script will start by gathering all the file necessary for the run, will copy them in the path "/home/ubuntu" and then execute the host.py script which will start the kernel.
At the end of the run it will remove all the files copied and upload the fan speed controller on the FPGA. The password for the root account will be asked three times.