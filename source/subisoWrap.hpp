#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP

#include "Parameters.hpp"
#include <hls_stream.h>
#include <ap_int.h>

void subisoWrapper(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        ap_uint<512> htb_buf[DDR_WIDTH],
        T_DDR res_buf[RES_WIDTH],

#ifdef DEBUG_INTERFACE
        unsigned int &debug_endpreprocess_s,
#endif

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        );

#endif
