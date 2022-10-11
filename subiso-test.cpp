#include <fstream>
#include <iostream>
#include <cstdio>
#include <unordered_map>

#include <hls_stream.h>
#include <ap_int.h>

#include "subgraphIsomorphism.hpp"

template<uint8_t V_ID_W, uint8_t V_L_W>
void stream_query(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end)
{

    /* Query data structure */
    uint8_t numQueryVertices = 0;       
    uint16_t numQueryEdges = 0;       

    /* Query files */
    std::ifstream fQueryOrd("data/queryOrder.txt");
    std::ifstream fQueryLab("data/queryLabels.txt");
    std::ifstream fQueryEdges("data/queryEdges.txt");
    std::string fLine{};

    std::unordered_map<unsigned long, unsigned long> vToLabel;

    /* Store labels */
    std::getline(fQueryLab, fLine);
    sscanf(fLine.c_str(), "%hhu", &numQueryVertices);	
    for(int count = 0; count < numQueryVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fQueryLab, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &node_t, &label_t);
        vToLabel.insert(std::make_pair(node_t, label_t));
    }

    /* Stream matching order */
    for(int count = 0; count < numQueryVertices; count++){ 
        unsigned long node_t;   
        getline (fQueryOrd, fLine);
        sscanf(fLine.c_str(), "%lu", &node_t);
        ap_uint<V_ID_W> node = node_t;

        stream_src.write(node);
        stream_end.write(false);
    }
    stream_end.write(true);

    /* Stream edges */
    std::getline(fQueryEdges, fLine);
    sscanf(fLine.c_str(), "%hhu", &numQueryEdges);	
    for(int count = 0; count < numQueryEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fQueryEdges, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<V_ID_W> nodesrc = nodesrc_t;
        ap_uint<V_ID_W> nodedst = nodedst_t;
        ap_uint<V_L_W> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<V_L_W> labeldst = vToLabel.at(nodedst_t);

        stream_src.write(nodesrc);
        stream_src_l.write(labelsrc);
        stream_dst.write(nodedst);
        stream_dst_l.write(labeldst);
        stream_end.write(false);
    }
    stream_end.write(true);

    fQueryLab.close();
    fQueryEdges.close();
    fQueryOrd.close();
}

template<uint8_t V_ID_W, uint8_t V_L_W>
void stream_datagraph(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end)
{


    uint32_t numDataVertices = 0;       
    uint32_t numDataEdges = 0;       

    /* Data graph files */
    std::ifstream fDataLab("data/dataLabels.txt");
    std::ifstream fDataEdges("data/dataEdges.txt");
    std::string fLine{};

    std::unordered_map<unsigned long, unsigned long> vToLabel;

    /* Store labels */
    std::getline(fDataLab, fLine);
    sscanf(fLine.c_str(), "%u", &numDataVertices);	
    for(int count = 0; count < numDataVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fDataLab, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &node_t, &label_t);
        vToLabel.insert(std::make_pair(node_t, label_t));
    }

    /* Stream edges */
    std::getline(fDataEdges, fLine);
    sscanf(fLine.c_str(), "%u", &numDataEdges);
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fDataEdges, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<V_ID_W> nodesrc = nodesrc_t;
        ap_uint<V_ID_W> nodedst = nodedst_t;
        ap_uint<V_L_W> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<V_L_W> labeldst = vToLabel.at(nodedst_t);

        stream_src.write(nodesrc);
        stream_src_l.write(labelsrc);
        stream_dst.write(nodedst);
        stream_dst_l.write(labeldst);
        stream_end.write(false);
    }
    stream_end.write(true);

    /* Rewind the file for 2nd stream */
    fDataEdges.clear();
    fDataEdges.seekg(0);

    /* Stream edges 2nd time */
    std::getline(fDataEdges, fLine);
    sscanf(fLine.c_str(), "%hhu", &numDataEdges);	
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fDataEdges, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<V_ID_W> nodesrc = nodesrc_t;
        ap_uint<V_ID_W> nodedst = nodedst_t;
        ap_uint<V_L_W> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<V_L_W> labeldst = vToLabel.at(nodedst_t);

        stream_src.write(nodesrc);
        stream_src_l.write(labelsrc);
        stream_dst.write(nodedst);
        stream_dst_l.write(labeldst);
        stream_end.write(false);
    }
    stream_end.write(true);

    fDataLab.close();
    fDataEdges.close();
}

int main()
{
    std::cout << "Start" << std::endl;

    hls::stream<ap_uint<VERTEX_WIDTH>> stream_src("src nodes");
    hls::stream<ap_uint<VERTEX_WIDTH>> stream_dst("dst nodes");
    hls::stream<ap_uint<VERTEX_WIDTH>> stream_out("result");
    hls::stream<ap_uint<LABEL_WIDTH>> stream_src_l("src labels");
    hls::stream<ap_uint<LABEL_WIDTH>> stream_dst_l("dst labels");
    hls::stream<bool> stream_end("bool ends");

    stream_query
        <VERTEX_WIDTH, 
        LABEL_WIDTH>
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             stream_end);

    stream_datagraph
        <VERTEX_WIDTH, 
        LABEL_WIDTH>
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             stream_end);
    
    subgraphIsomorphism
        <VERTEX_WIDTH, 
        LABEL_WIDTH,
        MAX_QUERY_VERTICES,
        MAX_TABLES,
        HASH_WIDTH,
        COUNTER_WIDTH>
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             stream_end,
             stream_out);

    int counter = 0;
    while(!stream_out.empty()){
        if (counter == 0){
            std::cout << "{ ";
        }
        std::cout << (int)stream_out.read() << " ";
        counter++;
        if (counter == 4){
            counter = 0;
            std::cout << "}" << std::endl;
        }
    }

    std::cout << "End" << std::endl;
    return 0;	
}
