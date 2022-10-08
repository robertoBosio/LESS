#include <fstream>
#include <ap_int.h>
#include <iostream>
#include <cstdio>
#include "subgraphIsomorphism.hpp"

void stream_query(
        hls::stream<ap_uint<V_ID_W>> stream_src,
        hls::stream<ap_uint<V_ID_W>> stream_dst,
        hls::stream<ap_uint<V_L_W>> stream_src_l,
        hls::stream<ap_uint<V_L_W>> stream_dst_l,
        hls::stream<ap_uint<bool>> stream_end)
{

	/* Query data structure */
	uint8_t numQueryVertices = 0;       
	uint16_t numQueryEdges = 0;       
    
    /* Query files */
	std::ifstream fQueryOrd("data/queryOrder.txt");
    std::ifstream fQueryLab("data/queryLables.txt");
	std::ifstream fQueryEdges("data/queryOffsets.txt");
	std::string fLine{};

    std::unordered_map<ap_uint<V_ID_W>, ap_uint<V_L_W>> vToLabel;

	/* Store labels */
	std::getline(fQueryLab, fLine);
	sscanf(fLine.c_str(), "%hhu", &numQueryVertices);	
	for(int count = 0; count < numQueryVertices; count++){    
		unsigned long node_t, label_t;
        std::getline(fQueryLab, fLine);
		sscanf(fLine.c_str(), "%lu %lu", &node_t, &label_t);
		ap_uint<V_ID_W> node = node_t;
        ap_uint<V_L_W> label = label_t;
		vToLabel.insert(std::make_pair(node, label));
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
        
        stream_src.write(nodesrc);
        stream_src_l.write(vToLabel.at(nodesrc));
        stream_dst.write(nodesrc);
        stream_dst_l.write(vToLabel.at(nodedst));
        stream_end.write(false);
	}
    stream_end.write(true);

    fQueryLab.close();
	fQueryEdges.close();
	fQueryOrd.close();
}

int main()
{
	std::cout << "Start" << std::endl;

    hls::stream<ap_uint<V_ID_W>> stream_src;
    hls::stream<ap_uint<V_ID_W>> stream_dst;
    hls::stream<ap_uint<V_L_W>> stream_src_l;
    hls::stream<ap_uint<V_L_W>> stream_dst_l;
    hls::stream<ap_uint<bool>> stream_end;

    streamQuery(
			stream_src,
			stream_dst,
			stream_src_l,
			stream_dst_l,
			stream_end);
    
	subgraphIsomorphism
        <VERTEX_WIDTH, 
        LABEL_WIDTH>
        (
			stream_src,
			stream_dst,
			stream_src_l,
			stream_dst_l,
			stream_end);

	std::cout << "End" << std::endl;
	return 0;	
}
