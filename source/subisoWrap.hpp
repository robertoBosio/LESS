#pragma once

#define HLS_STREAM_THREAD_SAFE
#define AP_INT_MAX_W 8192

#include "Parameters.hpp"
#include <ap_int.h>
#include <hls_stream.h>

#if SOFTWARE_PREPROC
#include "QueryVertex.hpp"
#include "Trie.hpp"

void subisoWrapper(
        row_t htb_buf[HASHTABLES_SPACE],
        bloom_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        unsigned long &diagnostic,
        QueryVertex qVertices0[MAX_QUERY_VERTICES], 
        QueryVertex qVertices1[MAX_QUERY_VERTICES],
        AdjHT hTables0[MAX_TABLES],
        AdjHT hTables1[MAX_TABLES],
#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */
#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */
        );

#else

void subisoWrapper(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf[HASHTABLES_SPACE],
        bloom_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned short numQueryEdges,
        const unsigned long numDataEdges,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        unsigned long &diagnostic,
#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */
#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */
        );

#endif /* SOFTWARE_PREPROC */
