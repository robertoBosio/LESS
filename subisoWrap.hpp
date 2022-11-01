#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP

#include "Parameters.hpp"
#include <hls_stream.h>
#include <ap_int.h>

void subisoWrap(
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_src,
        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_dst,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_src_l,
        hls::stream<ap_uint<LABEL_WIDTH>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<VERTEX_WIDTH_BIT>> &stream_out,
        hls::stream<bool> &stream_end_out);

#endif
