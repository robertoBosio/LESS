#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP

#include "Parameters.hpp"
#include <hls_stream.h>
#include <ap_int.h>

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

		ap_uint<VERTEX_WIDTH_BIT> res_buf[RES_WIDTH]);

#endif
