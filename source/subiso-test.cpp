#include "subgraphIsomorphism.hpp"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <unistd.h>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>

#include <ap_int.h>
#include <hls_stream.h>
#include "Parameters.hpp"
#include "debug.hpp"

#if SOFTWARE_PREPROC
#include "preprocess.hpp"
#endif /* SOFTWARE_PREPROC */

struct TestEntry {
    std::string querygraph;
    std::string golden;
    std::string h1;
    std::string h2;
};

template <size_t NODE_W,
         size_t LAB_W,
         size_t BURST_SIZE,
         size_t RESULT_SPACE,
         size_t MAX_QDATA>
void load_datagraphs(
        row_t *edge_buf,
        std::string datafile,
        unsigned long &dynfifo_space,
        unsigned long &numDataEdges)
{
    unsigned long numDataVertices;
    unsigned long edge_buf_p = 0;
    edge_t edge;
    
    /* Remove "../" to make paths correct */
    datafile = datafile.substr(3);
    
    std::ifstream fData(datafile);
    
    if (!fData.is_open()){
        std::cout << "Datagraph file opening failed.\n";
        return;
    }

    std::string fLine{};
    std::unordered_map<unsigned long, unsigned long> vToLabelData;

    std::getline(fData, fLine);
    sscanf(fLine.c_str(), "%*c %lu %lu", &numDataVertices, &numDataEdges);

    // Find space for the graph and align it to BURST_SIZE
    dynfifo_space = numDataEdges + MAX_QDATA;
    dynfifo_space = dynfifo_space - (dynfifo_space % BURST_SIZE) + BURST_SIZE;
    dynfifo_space = RESULT_SPACE - dynfifo_space;

    dynfifo_space = 100000;
    edge_buf_p = dynfifo_space;

    /* Store data labels */
    for(int count = 0; count < numDataVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fData, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu %*u", &node_t, &label_t);
        vToLabelData.insert(std::make_pair(node_t, label_t));
    }
    
    /* Stream edges */
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fData, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        edge.labelsrc = vToLabelData.at(nodesrc_t); 
        edge.labeldst = vToLabelData.at(nodedst_t);
        edge.src = nodesrc_t;
        edge.dst = nodedst_t;
        edge_buf[edge_buf_p++] = *((row_t*)&edge);
    }

    fData.close();
}

template <size_t NODE_W,
         size_t LAB_W,
         size_t MAX_QDATA>
void load_querygraphs(
        row_t *edge_buf,
        std::string queryfile,
        const unsigned long dynfifo_space,
        unsigned short &numQueryVertices,
        unsigned short &numQueryEdges,
        unsigned long numDataEdges)
{
    unsigned long numDataVertices;
    unsigned long edge_buf_p = numDataEdges + dynfifo_space;
    edge_t edge;
    
    /* Remove "../" to make paths correct */
    queryfile = queryfile.substr(3);
    
    /* Query file */
    std::ifstream fQuery(queryfile);

    // Get the current path
    char currentPath[FILENAME_MAX];
    if (getcwd(currentPath, sizeof(currentPath)) != nullptr) {
        std::cout << "Current Path: " << currentPath << std::endl;
    } else {
        std::cerr << "Error getting current path: " << strerror(errno) << std::endl;
        exit(-1);
    }
    std::cout << queryfile << std::endl; 
    if (!fQuery.is_open()){
        std::cout << "Query file opening failed.\n";
        exit(-1);
    }

    std::string fLine{};
    std::unordered_map<unsigned long, unsigned long> vToLabelQuery;

    /* Read vertices and edges cardinality */
    std::getline(fQuery, fLine);
    sscanf(fLine.c_str(), "%*c %hu %hu", &numQueryVertices, &numQueryEdges);

    assert(MAX_QDATA >= numQueryEdges + numQueryVertices);
    
    /* Store query labels */
    for (int count = 0; count < numQueryVertices; count++)
    {
        unsigned long node_t, label_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu %*u", &node_t, &label_t);
        vToLabelQuery.insert(std::make_pair(node_t, label_t));
    }

    /* Stream matching order 
     * RIGHT NOW NUMERICAL ORDERD */
    for(int count = 0; count < numQueryVertices; count++){ 
        edge.src = (ap_uint<NODE_W>)count;
        edge_buf[edge_buf_p++] = *((row_t*)&edge);
    }

    /* Stream edges */
    for(int count = 0; count < numQueryEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        edge.labelsrc = vToLabelQuery.at(nodesrc_t); 
        edge.labeldst = vToLabelQuery.at(nodedst_t);
        edge.src = nodesrc_t;
        edge.dst = nodedst_t;
        edge_buf[edge_buf_p++] = *((row_t*)&edge);
    }

    fQuery.close();
}

unsigned int subgraphIsomorphism_sw(){
    std::ifstream fGolden("data/golden.txt");
    unsigned int nEmbedd = 0;
    std::string fLine{};
    std::getline(fGolden, fLine);
    sscanf(fLine.c_str(), "%u", &nEmbedd);
    fGolden.close();
    return nEmbedd;
}

unsigned int countSol(
        int nQV,
        hls::stream<T_NODE> &stream_result)
{
    unsigned int nEmbedd = 0;
//    ap_uint<VERTEX_WIDTH_BIT> read;
    std::ofstream fres("data/resultkern.txt");

    T_NODE node;
    node = stream_result.read();

    while (!node.last){
    	for(int g = 0; g < nQV; g++){
            fres << node.data << " ";
            node = stream_result.read();
    	}
        fres << std::endl;
        nEmbedd++;
    }
    fres.close();
    return nEmbedd;	
}

int main()
{
#if COUNT_ONLY
    long unsigned int result;
#else
    hls::stream<T_NODE> result("results");
#endif


#if DEBUG_INTERFACE
    unsigned int debug_endpreprocess_s {0};
#endif

    // const unsigned char hash1_w = HASH_WIDTH_FIRST;
    // const unsigned char hash2_w = HASH_WIDTH_SECOND;
    // const unsigned char hash1_w = 10;
    // const unsigned char hash2_w = 6;
    // const unsigned long htb_size = (1UL << (hash1_w + hash2_w - (DDR_BIT - COUNTER_WIDTH)));
    unsigned short nQV = 0;
    unsigned short nQE = 0;
    unsigned long nDE = 0;
    unsigned int res_actual;
    unsigned int res_expected;
    unsigned long diagnostic;
    unsigned long dynfifo_space;
    bool flag = true;
    std::map<std::string, std::vector<TestEntry>> test;
    std::string prev_datagraph;
    unsigned long counters[25];
    unsigned int dynfifo_overflow;

    char cwd[100];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        std::cout << "Current working dir: " << cwd << std::endl;

    row_t *res_buf = (row_t*)malloc(RESULTS_SPACE * sizeof(row_t));
    
    if (!res_buf){
		std::cout << "Allocation failed." << std::endl;
		return -1;
	}

    std::cout << "Allocating " << 
        (unsigned long)((HASHTABLES_SPACE + RESULTS_SPACE) * sizeof(row_t) +
        BLOOM_SPACE * sizeof(bloom_t)) << " bytes." << std::endl;


    std::ifstream testfile("scripts/run_list.txt");
    if (!testfile) {
        std::cerr << "Error: Unable to open test.txt" << std::endl;
        return -1;
    }

    std::string line;
    while (std::getline(testfile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::string datagraph, querygraph, golden, h1, h2;
        std::istringstream iss(line);
        iss >> datagraph >> querygraph >> golden >> h1 >> h2;


        if (datagraph == prev_datagraph) {
            test[datagraph].push_back({querygraph, golden, h1, h2});
        } else {
            test[datagraph] = {{querygraph, golden, h1, h2}};
        }
        prev_datagraph = datagraph;
    }

    testfile.close();

    char datagraph_file[100], querygraph_file[100];
    row_t *htb_buf = (row_t *)calloc(HASHTABLES_SPACE, sizeof(row_t));
    if (!htb_buf)
    {
        std::cout << "Allocation failed." << std::endl;
        return -1;
    }

    bloom_t *bloom_p = (bloom_t *)calloc(BLOOM_SPACE, sizeof(bloom_t));
    if (!bloom_p)
    {
        std::cout << "Allocation failed." << std::endl;
        return -1;
    }

    for (const auto& entry : test) {
        const std::string& datagraph = entry.first;
        const std::vector<TestEntry>& entries = entry.second;

        std::cout << "Datagraph: " << datagraph << std::endl;

        // load graph
        load_datagraphs<
            VERTEX_WIDTH_BIT,
            LABEL_WIDTH,
            DYN_FIFO_BURST,
            RESULTS_SPACE,
            MAX_QUERYDATA>(
            res_buf,
            std::string(datagraph),
            dynfifo_space,
            nDE);

        for (const TestEntry &testEntry : entries)
        {

            // load query
            load_querygraphs<
                VERTEX_WIDTH_BIT,
                LABEL_WIDTH,
                MAX_QUERYDATA>(
                res_buf,
                std::string(testEntry.querygraph),
                dynfifo_space,
                nQV,
                nQE,
                nDE);

            std::cout << "  Querygraph: " << testEntry.querygraph << ", Golden: " << testEntry.golden
                      << ", H1: " << testEntry.h1 << ", H2: " << testEntry.h2 << std::endl;
            std::cout << "Words for dynamic fifo: " << dynfifo_space << std::endl;
            unsigned char h1 = stoi(testEntry.h1);
            unsigned char h2 = stoi(testEntry.h2);
            res_expected = stol(testEntry.golden);

#if SOFTWARE_PREPROC
            QueryVertex qVertices[MAX_QUERY_VERTICES];
            AdjHT hTables0[MAX_TABLES];
            AdjHT hTables1[MAX_TABLES];
            unsigned int n_candidate = 0;
            unsigned int start_candidate = 0;

            preprocess<row_t,
                       bloom_t,
                       EDGE_WIDTH,
                       COUNTER_WIDTH,
                       BLOOM_FILTER_WIDTH,
                       K_FUNCTIONS,
                       DDR_BIT,
                       VERTEX_WIDTH_BIT,
                       VERTEX_WIDTH,
                       HASH_LOOKUP3_BIT,
                       MAX_HASH_TABLE_BIT,
                       64,
                       LABEL_WIDTH,
                       DEFAULT_STREAM_DEPTH,
                       HASHTABLES_SPACE,
                       MAX_QUERY_VERTICES,
                       MAX_TABLES,
                       MAX_COLLISIONS>(res_buf,
                                       htb_buf,
                                       htb_buf,
                                       bloom_p,
                                       qVertices,
                                       hTables0,
                                       hTables1,
                                       dynfifo_space,
                                       n_candidate,
                                       start_candidate,
                                       nQV,
                                       nQE,
                                       nDE,
                                       h1,
                                       h2);

            subgraphIsomorphism(
                htb_buf,
                htb_buf,
                htb_buf,
                htb_buf,
                bloom_p,
                res_buf,
                nQV,
                h1,
                h2,
                dynfifo_space,
                dynfifo_overflow,
                n_candidate,
                start_candidate,
                qVertices,
                hTables0,
                hTables1,
#if DEBUG_INTERFACE
                debug_endpreprocess_s,
                counters[0],
                counters[1],
                counters[2],
                counters[3],
                counters[4],
                counters[5],
                counters[6],
                counters[7],
                counters[8],
                counters[9],
                counters[10],
                counters[11],
                counters[12],
                counters[13],
                counters[14],
                counters[15],
                counters[16],
                counters[17],
                counters[18],
                counters[19],
                counters[20],
                counters[21],
                counters[22],
                counters[23],
                counters[24],
#endif
                result);

#else
            subgraphIsomorphism(
                htb_buf,
                htb_buf,
                htb_buf,
                htb_buf,
                bloom_p,
                res_buf,
                nQV,
                nQE,
                nDE,
                h1,
                h2,
                dynfifo_space,
                dynfifo_overflow,
#if DEBUG_INTERFACE
                debug_endpreprocess_s,
                counters[0],
                counters[1],
                counters[2],
                counters[3],
                counters[4],
                counters[5],
                counters[6],
                counters[7],
                counters[8],
                counters[9],
                counters[10],
                counters[11],
                counters[12],
                counters[13],
                counters[14],
                counters[15],
                counters[16],
                counters[17],
                counters[18],
                counters[19],
                counters[20],
                counters[21],
                counters[22],
                counters[23],
                counters[24],
#endif
                result);
#endif /* SOFTWARE_PREPROC */

#if !COUNT_ONLY
            res_actual = countSol(
                nQV,
                result);
#else
            res_actual = result;
#endif

            std::cout << "Expected: " << res_expected << " actual: " << res_actual << std::endl;
            flag &= (res_actual == res_expected);
            for (int g = 0; g < HASHTABLES_SPACE; htb_buf[g++] = 0)
                ;
            for (int g = 0; g < BLOOM_SPACE; bloom_p[g++] = 0)
                ;
            for (int g = 0; g < 25; g++) {
                std::cout << counters[g] << std::endl;
            }
        }
    }
    free(htb_buf);
    free(bloom_p);
    free(res_buf);
    return !flag;
}
