#ifndef SUBISOWRAP_HPP
#define SUBISOWRAP_HPP
#define HLS_STREAM_THREAD_SAFE
#define AP_INT_MAX_W 8192

#include "Parameters.hpp"
#include <ap_int.h>
#include <hls_stream.h>
#include "hls_burst_maxi.h"
/* #include "hls_stream_sizeup.h" */

void subisoWrapper(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf[HASHTABLES_SPACE],
/* hls::burst_maxi<row_t> htb_buf1, */
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,

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
