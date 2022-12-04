#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

void subisoWrap(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        ap_uint<512> htb_buf0[DDR_WIDTH],
        ap_uint<512> htb_buf1[DDR_WIDTH],
        ap_uint<512> htb_buf2[DDR_WIDTH],
        ap_uint<512> htb_buf3[DDR_WIDTH],
        ap_uint<512> htb_buf4[DDR_WIDTH],

        hls::stream<T_NODE> &stream_out)
{

#pragma HLS INTERFACE m_axi port=htb_buf0 bundle=gmem0
#pragma HLS INTERFACE m_axi port=htb_buf1 bundle=gmem1
#pragma HLS INTERFACE m_axi port=htb_buf2 bundle=gmem2
#pragma HLS INTERFACE m_axi port=htb_buf3 bundle=gmem3
#pragma HLS INTERFACE m_axi port=htb_buf4 bundle=gmem4
#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3,htb_buf4 distance=0
#pragma HLS INTERFACE s_axilite port=return

#pragma HLS INTERFACE mode=axis port=stream_src
#pragma HLS INTERFACE mode=axis port=stream_dst
#pragma HLS INTERFACE mode=axis port=stream_src_l
#pragma HLS INTERFACE mode=axis port=stream_dst_l
#pragma HLS INTERFACE mode=axis port=stream_out

    subgraphIsomorphism
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             htb_buf0,
             htb_buf1,
             htb_buf2,
             htb_buf3,
             htb_buf4,
             stream_out);

}

