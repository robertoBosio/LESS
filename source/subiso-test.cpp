#include <fstream>
#include <iostream>
#include <cstdio>
#include <unordered_map>
#include <hls_stream.h>
#include <ap_int.h>

#include "subisoWrap.hpp"
#include "Parameters.hpp"

void stream_query(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        unsigned int &numQueryVertices)
{

    /* Query data structure */
    unsigned int numQueryEdges = 0;       

    /* Query files */
    std::ifstream fQueryOrd("data/queryOrder.txt");
    std::ifstream fQuery("data/querygraph.csv");
    std::string fLine{};

    std::unordered_map<unsigned long, unsigned long> vToLabel;

    /* Store labels */
    std::getline(fQuery, fLine);
    sscanf(fLine.c_str(), "%*c %u %u", &numQueryVertices, &numQueryEdges);	
    for(int count = 0; count < numQueryVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu %*u", &node_t, &label_t);
        vToLabel.insert(std::make_pair(node_t, label_t));
    }

    /* Stream matching order */
    for(int count = 0; count < numQueryVertices; count++){ 
        unsigned long node_t;   
        getline (fQueryOrd, fLine);
        sscanf(fLine.c_str(), "%lu", &node_t);
        ap_uint<VERTEX_WIDTH_BIT> node = node_t;

		T_NODE nodesrcif;
        nodesrcif.data = node;
        nodesrcif.last = (count == (numQueryVertices - 1));
        stream_src.write(nodesrcif);
    }

    /* Stream edges */
    for(int count = 0; count < numQueryEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fQuery, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

        T_LABEL labeldstif;
		T_LABEL labelsrcif;
		T_NODE nodedstif;
		T_NODE nodesrcif;
		labeldstif.data = labeldst;
		labelsrcif.data = labelsrc;
		nodedstif.data = nodedst;
		nodesrcif.data = nodesrc;
		nodesrcif.last = (count == (numQueryEdges - 1));

		stream_dst_l.write(labeldstif);
		stream_src_l.write(labelsrcif);
		stream_dst.write(nodedstif);
		stream_src.write(nodesrcif);
    }

    fQuery.close();
    fQueryOrd.close();
    //return numQueryVertices;
}

void stream_datagraph(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l)
{


    unsigned long numDataVertices = 0;       
    unsigned long numDataEdges = 0;       

    /* Data graph files */
    std::ifstream fData("data/datagraph.csv");
    std::string fLine{};

    std::unordered_map<unsigned long, unsigned long> vToLabel;

    /* Store labels */
    std::getline(fData, fLine);
    sscanf(fLine.c_str(), "%*c %lu %lu", &numDataVertices, &numDataEdges);	
    for(int count = 0; count < numDataVertices; count++){    
        unsigned long node_t, label_t;
        std::getline(fData, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu %*u", &node_t, &label_t);
        vToLabel.insert(std::make_pair(node_t, label_t));
    }

    /* Stream edges */
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fData, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

        T_LABEL labeldstif;
		T_LABEL labelsrcif;
		T_NODE nodedstif;
		T_NODE nodesrcif;
        labeldstif.data = labeldst;
		labelsrcif.data = labelsrc;
		nodedstif.data = nodedst;
		nodesrcif.data = nodesrc;
        nodesrcif.last = (count == (numDataEdges - 1));

		stream_dst_l.write(labeldstif);
		stream_src_l.write(labelsrcif);
		stream_dst.write(nodedstif);
		stream_src.write(nodesrcif);

    }

    /* Rewind the file for 2nd stream */
    fData.clear();
    fData.seekg(0);

    /* Skip vertices section */
    for(int count = 0; count < numDataVertices + 1; count++){    
        std::getline(fData, fLine);
    }

    /* Stream edges 2nd time */
    for(int count = 0; count < numDataEdges; count++){    
        unsigned long nodesrc_t, nodedst_t;
        std::getline(fData, fLine);
        sscanf(fLine.c_str(), "%*c %lu %lu", &nodesrc_t, &nodedst_t);
        ap_uint<VERTEX_WIDTH_BIT> nodesrc = nodesrc_t;
        ap_uint<VERTEX_WIDTH_BIT> nodedst = nodedst_t;
        ap_uint<LABEL_WIDTH> labelsrc = vToLabel.at(nodesrc_t);
        ap_uint<LABEL_WIDTH> labeldst = vToLabel.at(nodedst_t);

        T_LABEL labeldstif;
		T_LABEL labelsrcif;
		T_NODE nodedstif;
		T_NODE nodesrcif;
		labeldstif.data = labeldst;
		labelsrcif.data = labelsrc;
		nodedstif.data = nodedst;
		nodesrcif.data = nodesrc;
		nodesrcif.last = (count == (numDataEdges - 1));

		stream_dst_l.write(labeldstif);
		stream_src_l.write(labelsrcif);
		stream_dst.write(nodedstif);
		stream_src.write(nodesrcif);
    }

    fData.close();
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
    hls::stream<T_NODE> stream_src("src nodes");
    hls::stream<T_NODE> stream_dst("dst nodes");
    hls::stream<T_NODE> stream_out("result");
    hls::stream<T_LABEL> stream_src_l("src labels");
    hls::stream<T_LABEL> stream_dst_l("dst labels");

#ifdef COUNT_ONLY
    long unsigned int result;
#else
    hls::stream<T_NODE> result("results");
#endif


#ifdef DEBUG_INTERFACE
    unsigned int debug_endpreprocess_s {0};
#endif

    unsigned int nQV = 0;
    unsigned int res_actual;

    T_DDR *res_buf = (T_DDR*)malloc(RES_WIDTH * sizeof(T_DDR));
    
    if (!res_buf){
		std::cout << "Allocation failed." << std::endl;
		return -1;
	}


    std::cout << "Allocated " << 
        ((unsigned long)(DDR_WIDTH + RES_WIDTH) * DDR_WORD / 8 ) 
        << " bytes." << std::endl;
   
    for (int g = 0; g < 1; g++){ 
        
        ap_uint<DDR_WORD> *htb_buf = (ap_uint<DDR_WORD>*)calloc(DDR_WIDTH, sizeof(ap_uint<DDR_WORD>));
        if (!htb_buf){
            std::cout << "Allocation failed." << std::endl;
            return -1;
        }

        stream_query(
                stream_src,
                stream_dst,
                stream_src_l,
                stream_dst_l,
                nQV);

        stream_datagraph(
                stream_src,
                stream_dst,
                stream_src_l,
                stream_dst_l);

        subisoWrapper(
                stream_src,
                stream_dst,
                stream_src_l,
                stream_dst_l,
                htb_buf,
                res_buf,

#ifdef DEBUG_INTERFACE
                debug_endpreprocess_s,
#endif
                result);
       
#ifndef COUNT_ONLY 
        res_actual = countSol(
            nQV,
            result);
#else
        res_actual = result;
#endif

        free(htb_buf);
    }
    unsigned int res_expected = subgraphIsomorphism_sw();
    free(res_buf);
    return (res_actual != res_expected);
}
