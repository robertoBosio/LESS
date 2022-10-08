#include <ap_int.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"

#define DEBUG 1

template <uint8_t MAX_QV>
void buildTables(
        queryVertex qVer){


}

template <uint8_t V_ID_W, uint8_t V_L_W>
void subgraphIsomorphism(
        hls::stream<ap_uint<V_ID_W>> stream_src,
        hls::stream<ap_uint<V_ID_W>> stream_dst,
        hls::stream<ap_uint<V_L_W>> stream_src_l,
        hls::stream<ap_uint<V_L_W>> stream_dst_l,
        hls::stream<ap_uint<bool>> stream_end
        ){

#if DEBUG
    bool last = stream_end.read();
    while(!last){
        ap_uint<V_ID_W> node = stream_src.read();
        std::cout << (unsigned int)node << std::endl;
        last = stream_end.read();
    }

    last = stream_end.read();
    while(!last){
        ap_uint<V_ID_W> nodesrc = stream_src.read();
        ap_uint<V_ID_W> nodedst = stream_dst.read();
        ap_uint<V_L_W> labelsrc = stream_src_l.read();
        ap_uint<V_L_W> labeldst = stream_dst_l.read();
        std::cout << (unsigned int)nodesrc << "(" << (char)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (char)labeldst << ")" 
            << std::endl;
        last = stream_end.read();
    }
#endif

}
