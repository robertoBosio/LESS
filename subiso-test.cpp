#include <fstream>
#include <iostream>
#include <cstdio>
#include <unordered_map>

#include <hls_stream.h>
#include <ap_int.h>

#include "subisoWrap.hpp"
#include "Parameters.hpp"

int stream_query(
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_src,
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_dst,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_src_l,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_dst_l,
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
        ap_uint<VERTEX_WIDTH_BIT> node = node_t;

        stream_src.write(node);
        stream_end.write(false);
    }
    stream_end.write(true);

    /* Stream edges */
    std::getline(fQueryEdges, fLine);
    sscanf(fLine.c_str(), "%hu", &numQueryEdges);
    for(int count = 0; count < numQueryEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fQueryEdges, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

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
    return numQueryVertices;
}

void stream_datagraph(
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_src,
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_dst,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_src_l,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_dst_l,
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
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

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
    sscanf(fLine.c_str(), "%u", &numDataEdges);
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fDataEdges, fLine);
        sscanf(fLine.c_str(), "%lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

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
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_in,
        hls::stream<bool> &stream_end)
{
    bool last = stream_end.read();
    unsigned int nEmbedd = 0;
    ap_uint<VERTEX_WIDTH_BIT> temp;
    std::ofstream fres("data/resultkern.txt");
    while(!last){
        for (int g = 0; g < nQV; g++){
            temp = stream_in.read();
            fres << (int)temp << " "; 
        }  
        fres << std::endl;
        nEmbedd++;
        last = stream_end.read();
    }  
    fres.close();
    return nEmbedd;	
}

int main()
{
    std::cout << "Start" << std::endl;

    hls::stream<ap_uint<VERTEX_WIDTH_BIT>> stream_src("src nodes");
    hls::stream<ap_uint<VERTEX_WIDTH_BIT>> stream_dst("dst nodes");
    hls::stream<ap_uint<VERTEX_WIDTH_BIT>> stream_out("result");
    hls::stream<ap_uint<LABEL_WIDTH>> stream_src_l("src labels");
    hls::stream<ap_uint<LABEL_WIDTH>> stream_dst_l("dst labels");
    hls::stream<bool> stream_end_in("bool ends");
    hls::stream<bool> stream_end_out("end solutions");
    int nQV = 0;

    ap_uint<512> *htb_buf = (ap_uint<512>*)malloc(6000000);
    if (!htb_buf){
        std::cout << "Allocation failed." << std::endl;
        return -1;
    }

    nQV = stream_query(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end_in);

    stream_datagraph(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end_in);
    
    subisoWrap(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end_in,
            htb_buf,
            stream_out,
            stream_end_out);

    unsigned int res_actual = countSol(
            nQV,
            stream_out,
            stream_end_out);

    unsigned int res_expected = subgraphIsomorphism_sw();

    std::cout << "Expected: " << res_expected << ", Actual: " <<
        res_actual << std::endl;

    free(htb_buf);
    return (res_actual != res_expected);
}
