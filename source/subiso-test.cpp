#define HLS_STREAM_THREAD_SAFE
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <unordered_map>

#include <ap_int.h>
#include <hls_stream.h>
#include "subisoWrap.hpp"
#include "Parameters.hpp"
#include "debug.hpp"

#if SOFTWARE_PREPROC
#include "preprocess.hpp"
#endif /* SOFTWARE_PREPROC */

template <size_t NODE_W,
         size_t LAB_W,
         size_t GRAPH_SPACE>
edge_t *load_graphs(
        std::string queryfile,
        std::string datafile,
        unsigned short &numQueryVertices,
        unsigned short &numQueryEdges,
        unsigned long &numDataEdges)
{
    unsigned long numDataVertices;
    unsigned long edge_buf_p = 0;
    edge_t edge;

    /* Query file */
    std::ifstream fQuery(queryfile);
    std::ifstream fData(datafile);
    
    if (!fQuery.is_open() || !fData.is_open())
        return (edge_t*)nullptr;

    std::string fLine{};
    std::unordered_map<unsigned long, unsigned long> vToLabelQuery;
    std::unordered_map<unsigned long, unsigned long> vToLabelData;

    /* Read vertices and edges cardinality */
    std::getline(fQuery, fLine);
    sscanf(fLine.c_str(), "%*c %hu %hu", &numQueryVertices, &numQueryEdges);	
    std::getline(fData, fLine);
    sscanf(fLine.c_str(), "%*c %lu %lu", &numDataVertices, &numDataEdges);

    std::cout << numQueryVertices << " " << numQueryEdges << std::endl;
    std::cout << numDataVertices << " " << numDataEdges << std::endl;
    assert(GRAPH_SPACE > numDataEdges + numQueryVertices + numQueryEdges);

    edge_t *edge_buf = (edge_t*)malloc(GRAPH_SPACE * sizeof(edge_t));

    if (!edge_buf)
        return (edge_t*)nullptr;

    /* Store query labels */
    for(int count = 0; count < numQueryVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu %*u", &node_t, &label_t);
        vToLabelQuery.insert(std::make_pair(node_t, label_t));
    }
    
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
        edge.src = (nodesrc_t << LAB_W) | vToLabelData.at(nodesrc_t);
        edge.dst = (nodedst_t << LAB_W) | vToLabelData.at(nodedst_t);
        edge_buf[edge_buf_p++] = edge;
    }

    /* Stream matching order 
     * RIGHT NOW NUMERICAL ORDERD */
    for(int count = 0; count < numQueryVertices; count++){ 
        edge.src = (ap_uint<NODE_W>)count;
        edge_buf[edge_buf_p++] = edge;
    }

    /* Stream edges */
    for(int count = 0; count < numQueryEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        edge.src = (nodesrc_t << LAB_W) | vToLabelQuery.at(nodesrc_t);
        edge.dst = (nodedst_t << LAB_W) | vToLabelQuery.at(nodedst_t);
        edge_buf[edge_buf_p++] = edge;
    }

    fData.close();
    fQuery.close();
    return edge_buf;
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

    unsigned short nQV = 0;
    unsigned short nQE = 0;
    unsigned long nDE = 0;
    unsigned int res_actual;
    unsigned int res_expected;
    unsigned long diagnostic;
    bool flag = true;


    row_t *res_buf = (row_t*)malloc(RESULTS_SPACE * sizeof(row_t));
    
    if (!res_buf){
		std::cout << "Allocation failed." << std::endl;
		return -1;
	}


    std::cout << "Allocating " << 
        (unsigned long)((HASHTABLES_SPACE + RESULTS_SPACE) * sizeof(row_t) +
        BLOOM_SPACE * sizeof(bloom_t) +
        GRAPHS_SPACE * sizeof(edge_t)) << " bytes." << std::endl;
 
    std::ifstream fTest("data/test.csv");
    std::string fLine{};
    char datagraph_file[100], querygraph_file[100];
    std::getline(fTest, fLine);
    
    while (!fTest.eof()){
        if (fLine.c_str()[0] == '#'){
            std::getline(fTest, fLine);
            continue;
        }
        
        sscanf(fLine.c_str(), "%s %s %u", datagraph_file, querygraph_file, &res_expected);	
        row_t *htb_buf = (row_t*)calloc(HASHTABLES_SPACE, sizeof(row_t));
        if (!htb_buf){
            std::cout << "Allocation failed." << std::endl;
            return -1;
        }

        bloom_t *bloom_p = (bloom_t*)calloc(BLOOM_SPACE, sizeof(bloom_t));
        if (!bloom_p){
            std::cout << "Allocation failed." << std::endl;
            return -1;
        }

        edge_t *edge_buf = load_graphs<
            VERTEX_WIDTH_BIT,
            LABEL_WIDTH,
            GRAPHS_SPACE>(
                std::string(querygraph_file),
                std::string(datagraph_file),
                nQV,
                nQE,
                nDE);

        if (!edge_buf)
            return -1;

#if SOFTWARE_PREPROC 
        QueryVertex qVertices0[MAX_QUERY_VERTICES]; 
        QueryVertex qVertices1[MAX_QUERY_VERTICES];
        AdjHT hTables0[MAX_TABLES];
        AdjHT hTables1[MAX_TABLES];
        
        preprocess<row_t,
            bloom_t,
            EDGE_WIDTH,
            COUNTER_WIDTH,
            BLOOM_FILTER_WIDTH,
            K_FUNCTIONS,
            DDR_BIT,
            VERTEX_WIDTH_BIT,
            LABEL_WIDTH,
            HASH_WIDTH_FIRST,
            HASH_WIDTH_SECOND,
            DEFAULT_STREAM_DEPTH,
            HASHTABLES_SPACE,
            MAX_QUERY_VERTICES,
            MAX_TABLES>(
                    edge_buf,
                    htb_buf,
                    bloom_p,
                    qVertices0,
                    qVertices1,
                    hTables0,
                    hTables1,
                    nQV,
                    nQE,
                    nDE);
        
        subisoWrapper(
                htb_buf,
                bloom_p,
                res_buf,
                nQV,
                diagnostic,
                qVertices0,
                qVertices1,
                hTables0,
                hTables1,
#if DEBUG_INTERFACE
                debug_endpreprocess_s,
#endif
                result);

#else
        subisoWrapper(
                edge_buf,
                htb_buf,
                bloom_p,
                res_buf,
                nQV,
                nQE,
                nDE,
                diagnostic,
#if DEBUG_INTERFACE
                debug_endpreprocess_s,
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

        free(htb_buf);
        free(edge_buf);
        free(bloom_p);
        std::cout << "Expected: " << res_expected << 
            " actual: " << res_actual << std::endl;
        flag &= (res_actual == res_expected);
        std::getline(fTest, fLine);
    }

    free(res_buf);
    return !flag;
}
