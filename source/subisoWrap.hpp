#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP

#include "Parameters.hpp"
#ifndef __SYNTHESIS__
#include "myhls_stream.h"
#else
#include <hls_stream.h>
#endif
#include <ap_int.h>

void subisoWrapper(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        ap_uint<512> htb_buf[DDR_WIDTH],

		hls::stream<T_NODE> &stream_result);

#endif
