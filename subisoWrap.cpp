#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

void subisoWrap(
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_src,
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_dst,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_src_l,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_dst_l,
        hls::stream<bool> &stream_end_in,
        ap_uint<512> htb_buf[DDR_WIDTH],

        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_out,
        hls::stream<bool> &stream_end_out){

#pragma HLS INTERFACE m_axi port=htb_buf bundle=gmem0
//#pragma HLS INTERFACE m_axi port=htb_buf bundle=gmem1
//#pragma HLS INTERFACE m_axi port=htb_buf2 bundle=gmem0_2
//#pragma HLS INTERFACE m_axi port=htb_buf3 bundle=gmem0_3
//#pragma HLS INTERFACE m_axi port=htb_buf4 bundle=gmem0_4

/*
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_src
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_dst
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_src_l
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_dst_l
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_end_in
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_out
#pragma HLS INTERFACE mode=ap_fifo depth=1 port=stream_end_out
*/

    subgraphIsomorphism
        <VERTEX_WIDTH_BIT, 
        LABEL_WIDTH,
        MAX_QUERY_VERTICES,
        MAX_TABLES,
        HASH_WIDTH_FIRST,
        HASH_WIDTH_SECOND,
        COUNTER_WIDTH,
        EDGE_WIDTH,
        MAX_COLLISIONS>
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             stream_end_in,
             htb_buf,
             stream_out,
             stream_end_out);

}

