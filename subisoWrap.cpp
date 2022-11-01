#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

void subisoWrap(
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_src,
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_dst,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_src_l,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_dst_l,
        hls::stream<bool> &stream_end_in,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_out,
        hls::stream<bool> &stream_end_out){


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
