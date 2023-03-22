#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP
#define HLS_STREAM_THREAD_SAFE

#include "Parameters.hpp"
#include <ap_int.h>
#include <hls_stream.h>
#include "hls_burst_maxi.h"
/* #include "hls_stream_sizeup.h" */

void subisoWrapper(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        T_DDR htb_buf[DDR_WIDTH],
/* hls::burst_maxi<T_DDR> htb_buf1, */
        T_DDR res_buf[RES_WIDTH],

#ifdef DEBUG_INTERFACE
        unsigned int &debif_endpreprocess,
#endif

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        );

#endif
